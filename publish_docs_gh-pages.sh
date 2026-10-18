#!/usr/bin/env bash
# ======================================================================================
# Publish a pre-built static site (HTML/plots) to GitHub Pages without committing
# build artifacts to your main repo history.
#
# OVERVIEW
# --------
# - You run your heavy pipeline locally/cluster and it produces web-ready files in a
#   directory (default: ./results/publish_docs relative to your *current working directory*).
# - This script publishes that directory to a dedicated Pages branch as a single
#   force-pushed snapshot commit with no parent, so the Pages branch ALWAYS contains
#   exactly one commit and your main branch and history stay clean.
# - The directory is either a finished site with a top-level index.html, or a VitePress
#   homepage source tree with a top-level index.md that the workflow in its own
#   .github/workflows builds whenever the Pages branch is pushed.
# - Configure GitHub Pages to serve from: Settings -> Pages -> "Deploy from a branch"
#   -> Branch: (the Pages branch, default: gh-pages), Folder: /.
#
# HOW THE SNAPSHOT IS BUILT
# -------------------------
# The site is read where it already sits; it is never copied. `git` is pointed at the
# site directory as its work tree, lists its files, and stages them into a temporary
# index; the resulting tree is committed and pushed directly:
#
#   git --work-tree=<site> ls-files --others --exclude-per-directory=.gitignore \
#     | git update-index --add --stdin           (into GIT_INDEX_FILE=<temp index>)
#   git write-tree                               -> tree of exactly the site's files
#   git commit-tree <tree>                       -> commit with no parent
#   git push -f <remote> <commit>:<pages branch>
#
# Consequences worth knowing:
# - The only temporary artifact is that index (tens of KB), written inside the git
#   directory so that the repository's own filesystem is the only one needing free
#   space. A site of any size needs no scratch space, and TMPDIR is not used.
# - The only ignore rules applied are the site's own .gitignore files, which is how a
#   VitePress homepage keeps node_modules and its build caches out of the Pages branch.
#   With the work tree redirected to the site directory the repository's top-level
#   .gitignore is not consulted, and $GIT_DIR/info/exclude and core.excludesFile are
#   deliberately not read, so nothing in your own git setup can silently drop files out
#   of the published site.
# - The tree is built from scratch every run, so the snapshot is exact: a file deleted
#   from the site directory disappears from the published branch.
# - No local branch is created or deleted, and no worktree is added; the commit is
#   pushed by object name. Your refs are left exactly as they were.
# - The blobs are written into the repository's object store, as any commit's are. Each
#   publish supersedes the last, so the previous snapshot's objects become unreachable
#   and accumulate locally; `git gc --prune=now` reclaims them. They are never pushed:
#   `git push` sends only objects reachable from the ref being pushed.
#
# SAFETY / STRICT MODE
# --------------------
# - We use `set -Eeuo pipefail`:
#     - -E: traps propagate through functions/subshells
#     - -e: exit on any non-zero command
#     - -u: treat unset vars as errors
#     - -o pipefail: pipelines fail if *any* command fails
# - We narrow IFS and install traps so the temporary index is removed on exit. SIGKILL
#   cannot be trapped; a leftover index is tens of KB inside the git directory, and a
#   stale one is swept at the start of the next run.
# - The push is a force-push, so the script refuses to run unless the staged snapshot
#   has a top-level index.html or index.md, rather than replacing a live site with
#   nothing.
#
# USAGE
# -----
#   ./publish_docs_gh-pages.sh
#
# Examples:
#   PUBLISH_DOCS_GH_PAGES_SITE_DIR=out/site \
#   PUBLISH_DOCS_GH_PAGES_BRANCH=gh-pages \
#   PUBLISH_DOCS_GH_PAGES_REMOTE=origin \
#     ./publish_docs_gh-pages.sh
#
# ENV / ARGS
# ----------
# - PUBLISH_DOCS_GH_PAGES_SITE_DIR  (default: results/publish_docs; **relative to caller's CWD**)
#   Must contain a top-level index.html or index.md
# - PUBLISH_DOCS_GH_PAGES_BRANCH    (default: gh-pages)
# - PUBLISH_DOCS_GH_PAGES_REMOTE    (default: origin)
#
# REQUIREMENTS
# ------------
# - Run from *inside* your Git repo (any subdirectory is fine).
# - `git` >= 2.13 (for `rev-parse --absolute-git-dir`). No other tools are needed.
# - Push access to the remote.
# - Keep within GitHub Pages practical limits (site <= ~1 GB; single file <100 MB).
# ======================================================================================

set -Eeuo pipefail
IFS=$'\n\t'

# ------------------------- Configuration (env-overridable) -------------------------
PUBLISH_DOCS_GH_PAGES_SITE_DIR="${PUBLISH_DOCS_GH_PAGES_SITE_DIR:-results/publish_docs}"  # relative to caller CWD
PUBLISH_DOCS_GH_PAGES_BRANCH="${PUBLISH_DOCS_GH_PAGES_BRANCH:-gh-pages}"
PUBLISH_DOCS_GH_PAGES_REMOTE="${PUBLISH_DOCS_GH_PAGES_REMOTE:-origin}"

# ------------------------- Helpers -------------------------
_abort() { echo "ERROR: $*" >&2; exit 1; }
_info()  { printf '%s\n' "$*"; }

# Resolve site dir relative to caller's CWD -> absolute path
CALLER_CWD="$(pwd -P)"
if [[ "${PUBLISH_DOCS_GH_PAGES_SITE_DIR}" = /* ]]; then
  SITE_DIR_ABS="${PUBLISH_DOCS_GH_PAGES_SITE_DIR}"
else
  SITE_DIR_ABS="${CALLER_CWD}/${PUBLISH_DOCS_GH_PAGES_SITE_DIR}"
fi

# ------------------------- Step 1: Repo checks -------------------------
_info "[1/8] Verifying Git repository context…"
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  _abort "Run this from inside your Git repository."
fi
REPO_ROOT="$(git rev-parse --show-toplevel)"
# The git directory is not always "$REPO_ROOT/.git": in a submodule or a linked worktree
# that path is a file pointing elsewhere.
GIT_DIR_ABS="$(git rev-parse --absolute-git-dir)"
cd "$REPO_ROOT"

# ------------------------- Step 2: Git config checks -------------------------
_info "[2/8] Verifying Git user configuration…"
if ! git config user.name >/dev/null 2>&1; then
  _abort "Git user.name is not configured. Run: git config user.name 'Your Name'"
fi
if ! git config user.email >/dev/null 2>&1; then
  _abort "Git user.email is not configured. Run: git config user.email 'your@email.com'"
fi

# ------------------------- Step 3: Remote checks -------------------------
_info "[3/8] Checking remote '${PUBLISH_DOCS_GH_PAGES_REMOTE}'…"
if ! git remote get-url "${PUBLISH_DOCS_GH_PAGES_REMOTE}" >/dev/null 2>&1; then
  _abort "Remote '${PUBLISH_DOCS_GH_PAGES_REMOTE}' not found. Add with: git remote add ${PUBLISH_DOCS_GH_PAGES_REMOTE} <url>"
fi

# ------------------------- Step 4: Default branch guard -------------------------
_info "[4/8] Protecting the remote's default branch…"
DEFAULT_BRANCH="$(
  git ls-remote --symref "${PUBLISH_DOCS_GH_PAGES_REMOTE}" HEAD 2>/dev/null \
    | awk '$1 == "ref:" && $3 == "HEAD" { sub(/^refs\/heads\//, "", $2); print $2; exit }'
)"
if [[ -n "$DEFAULT_BRANCH" && "${PUBLISH_DOCS_GH_PAGES_BRANCH}" == "$DEFAULT_BRANCH" ]]; then
  _abort "Refusing to publish to the remote's default branch '${DEFAULT_BRANCH}'. Set PUBLISH_DOCS_GH_PAGES_BRANCH=gh-pages (or another non-default branch)."
fi
case "${PUBLISH_DOCS_GH_PAGES_BRANCH}" in
  main|master|develop)
    _abort "Refusing to publish to '${PUBLISH_DOCS_GH_PAGES_BRANCH}'. Choose a dedicated Pages branch (e.g., 'gh-pages')."
  ;;
esac

# ------------------------- Step 5: Source directory checks -------------------------
_info "[5/8] Validating site directory at '${SITE_DIR_ABS}'…"
if [[ ! -d "$SITE_DIR_ABS" ]]; then
  _abort "Directory not found: ${SITE_DIR_ABS}. Did your pipeline write results there?"
fi
if [[ ! -f "$SITE_DIR_ABS/index.html" && ! -f "$SITE_DIR_ABS/index.md" ]]; then
  _abort "Missing index.html or index.md at the top level of ${SITE_DIR_ABS}."
fi

# ------------------------- Step 6: Stage the site into a temporary index -------------
_info "[6/8] Staging '${SITE_DIR_ABS}' into a temporary index…"

# Sweep scratch left behind by a run that was killed before its trap could fire. Only
# entries over a day old, so a concurrent publish is never disturbed.
find "$GIT_DIR_ABS" -maxdepth 1 -name 'publish-pages.*' -mtime +0 -exec rm -rf {} + 2>/dev/null || true

# `git` will not read a zero-byte index, and `mktemp` has to create what it reserves, so
# reserve a directory and let `git` create the index inside it.
TMP_DIR="$(mktemp -d "${GIT_DIR_ABS}/publish-pages.XXXXXX")"
TMP_INDEX="${TMP_DIR}/index"
cleanup() {
  [[ -n "${TMP_DIR:-}" ]] && rm -rf "$TMP_DIR"
  return 0
}
trap cleanup EXIT
trap 'echo "Publish failed on line $LINENO" >&2; exit 1' ERR
trap 'exit 130' INT
trap 'exit 143' TERM HUP

# List every file in the site that its own .gitignore files do not exclude. Without
# `--exclude-standard`, `$GIT_DIR/info/exclude` and `core.excludesFile` are not read.
# Both commands use the temporary index, so the listing is against an empty index rather
# than the repository's own, and the paths are resolved against the site, hence the `-C`.
GIT_INDEX_FILE="$TMP_INDEX" \
  git -C "$SITE_DIR_ABS" --git-dir="$GIT_DIR_ABS" --work-tree="$SITE_DIR_ABS" \
    ls-files -z --others --exclude-per-directory=.gitignore \
  | GIT_INDEX_FILE="$TMP_INDEX" \
    git -C "$SITE_DIR_ABS" --git-dir="$GIT_DIR_ABS" --work-tree="$SITE_DIR_ABS" \
      update-index --add -z --stdin

# Allow files with leading underscores, etc. Added to the index rather than written into
# the site directory, which belongs to the pipeline that generated it.
EMPTY_BLOB="$(git --git-dir="$GIT_DIR_ABS" hash-object -w --stdin </dev/null)"
GIT_INDEX_FILE="$TMP_INDEX" \
  git --git-dir="$GIT_DIR_ABS" update-index --add --cacheinfo "100644,${EMPTY_BLOB},.nojekyll"

TREE="$(GIT_INDEX_FILE="$TMP_INDEX" git --git-dir="$GIT_DIR_ABS" write-tree)"

# ------------------------- Step 7: Build the snapshot commit -------------------------
_info "[7/8] Building snapshot commit…"

# The push below is forced, so verify there is a site in the tree before replacing the
# published one with it.
if ! git --git-dir="$GIT_DIR_ABS" cat-file -e "${TREE}:index.html" 2>/dev/null \
  && ! git --git-dir="$GIT_DIR_ABS" cat-file -e "${TREE}:index.md" 2>/dev/null; then
  _abort "Refusing to publish: the staged snapshot has no top-level index.html or index.md, so force-pushing it would destroy the published site."
fi
N_FILES="$(git --git-dir="$GIT_DIR_ABS" ls-tree -r --name-only "$TREE" | wc -l | tr -d '[:space:]')"
_info "• Snapshot holds ${N_FILES} files."

SRC_COMMIT="$(git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
DATE_UTC="$(date -u +'%Y-%m-%d %H:%M:%S UTC')"
COMMIT="$(
  git --git-dir="$GIT_DIR_ABS" commit-tree "$TREE" \
    -m "Publish site snapshot (${DATE_UTC}) from ${SRC_COMMIT} [source: ${SITE_DIR_ABS}]"
)"

# ------------------------- Step 8: Push -------------------------
_info "[8/8] Force-pushing to '${PUBLISH_DOCS_GH_PAGES_REMOTE}/${PUBLISH_DOCS_GH_PAGES_BRANCH}'…"
# The commit is pushed by object name, so no ref is created and this run touches nothing
# outside its own scratch directory. That is what lets two publishes overlap without
# interfering; `git gc` cannot collect the commit from under us either, as it only prunes
# unreferenced objects older than `gc.pruneExpire`, two weeks by default.
git push -f "${PUBLISH_DOCS_GH_PAGES_REMOTE}" "${COMMIT}:refs/heads/${PUBLISH_DOCS_GH_PAGES_BRANCH}"

_info "✓ Published '${SITE_DIR_ABS}' to ${PUBLISH_DOCS_GH_PAGES_REMOTE}/${PUBLISH_DOCS_GH_PAGES_BRANCH} as ${COMMIT}."
_info "      If not already configured, set GitHub Pages to serve from that branch (Settings → Pages)."

# End of script
