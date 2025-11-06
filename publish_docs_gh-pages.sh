#!/usr/bin/env bash
# ======================================================================================
# Publish a pre-built static site (HTML/plots) to GitHub Pages without committing
# build artifacts to your main repo history.
#
# OVERVIEW
# --------
# - You run your heavy pipeline locally/cluster and it produces web-ready files in a
#   directory (default: ./results/publish_docs relative to your *current working directory*).
# - This script publishes that directory to a dedicated Pages branch via a temporary
#   Git worktree and a single force-pushed snapshot commit.
# - Each publish *completely replaces* the branch with one fresh orphan commit that has
#   no parent history. The gh-pages branch will ALWAYS contain exactly one commit.
# - Your main branch and repo history stay clean, and the Pages branch never accumulates
#   history from previous publishes.
# - Configure GitHub Pages to serve from: Settings → Pages → "Deploy from a branch"
#   → Branch: (the Pages branch, default: gh-pages), Folder: /.
#
# WHAT THE KEY COMMANDS DO
# ------------------------
# - `git fetch --prune <remote>`:
#     Refreshes local knowledge of the remote’s refs and prunes stale ones so we can
#     accurately detect whether the Pages branch exists remotely.
#
# - Temporary worktree (TMP_DIR) via `git worktree`:
#     Creates an *isolated* working directory attached to the Pages branch. This keeps
#     your primary working tree untouched and gives us a clean target to mirror files.
#
# - `rsync -a --delete`:
#     Copies the site contents into the worktree and *deletes* anything that’s no longer
#     present in the source. The deployed site becomes an exact snapshot of your outputs.
#
# SAFETY / STRICT MODE
# --------------------
# - We use `set -Eeuo pipefail`:
#     - -E: traps propagate through functions/subshells
#     - -e: exit on any non-zero command
#     - -u: treat unset vars as errors
#     - -o pipefail: pipelines fail if *any* command fails
# - We narrow IFS and install traps so the temporary worktree is cleaned up on exit.
#   (SIGKILL cannot be trapped, but that’s rare; `git worktree prune` handles stragglers.)
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
# - PUBLISH_DOCS_GH_PAGES_SITE_DIR  (default: results/publish_docs; **relative to caller’s CWD**)
#   Must contain a top-level index.html
# - PUBLISH_DOCS_GH_PAGES_BRANCH    (default: gh-pages)
# - PUBLISH_DOCS_GH_PAGES_REMOTE    (default: origin)
#
# REQUIREMENTS
# ------------
# - Run from *inside* your Git repo (any subdirectory is fine).
# - Push access to the remote.
# - Keep within GitHub Pages practical limits (site ≤ ~1 GB; single file <100 MB).
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

# Resolve site dir relative to caller’s CWD → absolute path
CALLER_CWD="$(pwd -P)"
if [[ "${PUBLISH_DOCS_GH_PAGES_SITE_DIR}" = /* ]]; then
  SITE_DIR_ABS="${PUBLISH_DOCS_GH_PAGES_SITE_DIR}"
else
  SITE_DIR_ABS="${CALLER_CWD}/${PUBLISH_DOCS_GH_PAGES_SITE_DIR}"
fi

# ------------------------- Step 1: Repo checks -------------------------
_info "[1/10] Verifying Git repository context…"
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  _abort "Run this from inside your Git repository."
fi
REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

# ------------------------- Step 2: Git config checks -------------------------
_info "[2/10] Verifying Git user configuration…"
if ! git config user.name >/dev/null 2>&1; then
  _abort "Git user.name is not configured. Run: git config user.name 'Your Name'"
fi
if ! git config user.email >/dev/null 2>&1; then
  _abort "Git user.email is not configured. Run: git config user.email 'your@email.com'"
fi

# ------------------------- Step 3: Remote checks -------------------------
_info "[3/10] Checking remote '${PUBLISH_DOCS_GH_PAGES_REMOTE}'…"
if ! git remote get-url "${PUBLISH_DOCS_GH_PAGES_REMOTE}" >/dev/null 2>&1; then
  _abort "Remote '${PUBLISH_DOCS_GH_PAGES_REMOTE}' not found. Add with: git remote add ${PUBLISH_DOCS_GH_PAGES_REMOTE} <url>"
fi

# ------------------------- Step 4: Default branch guard -------------------------
_info "[4/10] Protecting the remote's default branch…"
DEFAULT_BRANCH="$(git remote show "${PUBLISH_DOCS_GH_PAGES_REMOTE}" | sed -n 's/^\s*HEAD branch: //p')"
if [[ -n "$DEFAULT_BRANCH" && "${PUBLISH_DOCS_GH_PAGES_BRANCH}" == "$DEFAULT_BRANCH" ]]; then
  _abort "Refusing to publish to the remote's default branch '${DEFAULT_BRANCH}'. Set PUBLISH_DOCS_GH_PAGES_BRANCH=gh-pages (or another non-default branch)."
fi
case "${PUBLISH_DOCS_GH_PAGES_BRANCH}" in
  main|master|develop)
    _abort "Refusing to publish to '${PUBLISH_DOCS_GH_PAGES_BRANCH}'. Choose a dedicated Pages branch (e.g., 'gh-pages')."
  ;;
esac

# ------------------------- Step 5: Source directory checks -------------------------
_info "[5/10] Validating site directory at '${SITE_DIR_ABS}'…"
if [[ ! -d "$SITE_DIR_ABS" ]]; then
  _abort "Directory not found: ${SITE_DIR_ABS}. Did your pipeline write results there?"
fi
if [[ ! -f "$SITE_DIR_ABS/index.html" && ! -f "$SITE_DIR_ABS/index.md" ]]; then
  _abort "Missing index.html or index.md at the top level of ${SITE_DIR_ABS}."
fi

# ------------------------- Step 6: Refresh remote refs -------------------------
_info "[6/10] Fetching remote refs (prune) from '${PUBLISH_DOCS_GH_PAGES_REMOTE}'…"
git fetch --prune "${PUBLISH_DOCS_GH_PAGES_REMOTE}"

# Clean up any stale worktrees and delete local branch to ensure clean orphan creation
_info "• Cleaning up any existing worktrees for '${PUBLISH_DOCS_GH_PAGES_BRANCH}'…"
git worktree prune
# Remove any worktrees still associated with the Pages branch
while IFS= read -r worktree_path; do
  if [[ -n "$worktree_path" ]]; then
    _info "  - Removing stale worktree at: $worktree_path"
    git worktree remove --force "$worktree_path" 2>/dev/null || true
  fi
done < <(git worktree list --porcelain | grep -A 2 "branch refs/heads/${PUBLISH_DOCS_GH_PAGES_BRANCH}" | grep "^worktree " | cut -d' ' -f2-)

# Now safe to delete the local branch if it exists
if git show-ref --verify --quiet "refs/heads/${PUBLISH_DOCS_GH_PAGES_BRANCH}"; then
  _info "• Deleting existing local branch '${PUBLISH_DOCS_GH_PAGES_BRANCH}' to create fresh orphan…"
  git branch -D "${PUBLISH_DOCS_GH_PAGES_BRANCH}"
fi

# ------------------------- Step 7: Create temporary worktree -------------------------
_info "[7/10] Creating temporary worktree for '${PUBLISH_DOCS_GH_PAGES_BRANCH}'…"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/publish-pages.XXXXXX")"
cleanup() {
  if git -C "$REPO_ROOT" worktree list | grep -Fq " $TMP_DIR "; then
    git -C "$REPO_ROOT" worktree remove --force "$TMP_DIR" || true
  fi
  [[ -d "$TMP_DIR" ]] && rm -rf "$TMP_DIR" || true
}
trap cleanup EXIT
trap 'echo "Publish failed on line $LINENO" >&2; exit 1' ERR
trap 'exit 130' INT
trap 'exit 143' TERM HUP

# Always create a fresh orphan branch (no history)
_info "• Creating fresh orphan branch (no history)…"
git worktree add --no-checkout "$TMP_DIR"
(
  cd "$TMP_DIR"
  # Create orphan branch with no parent commits
  git checkout --orphan "${PUBLISH_DOCS_GH_PAGES_BRANCH}"
  # Ensure index starts empty (working tree will be populated by rsync)
  git rm -r --cached . >/dev/null 2>&1 || true
)

# ------------------------- Step 8: Rsync snapshot -------------------------
_info "[8/10] Syncing '${SITE_DIR_ABS}' → worktree (exact snapshot)…"
rsync -a --delete --exclude='.git' "${SITE_DIR_ABS}/" "${TMP_DIR}/"
touch "${TMP_DIR}/.nojekyll"   # allow files with leading underscores, etc.

# ------------------------- Step 9: Commit & push -------------------------
_info "[9/10] Committing and force-pushing to '${PUBLISH_DOCS_GH_PAGES_REMOTE}/${PUBLISH_DOCS_GH_PAGES_BRANCH}'…"
(
  cd "$TMP_DIR"
  git add -A

  SRC_COMMIT="$(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
  DATE_UTC="$(date -u +'%Y-%m-%d %H:%M:%S UTC')"

  # Create a single snapshot commit (orphan branch has no history)
  if ! git diff --cached --quiet; then
    git commit --quiet -m "Publish site snapshot (${DATE_UTC}) from ${SRC_COMMIT} [source: ${SITE_DIR_ABS}]"
  else
    _info "• Nothing to publish; creating empty commit so branch exists…"
    git commit --quiet --allow-empty -m "Publish empty snapshot (${DATE_UTC})"
  fi

  # Force-push to completely replace the remote branch (no history preserved)
  git push -f -u "${PUBLISH_DOCS_GH_PAGES_REMOTE}" "HEAD:${PUBLISH_DOCS_GH_PAGES_BRANCH}"
)

# ------------------------- Step 10: Done -------------------------
_info "[10/10] ✓ Published '${SITE_DIR_ABS}' to ${PUBLISH_DOCS_GH_PAGES_REMOTE}/${PUBLISH_DOCS_GH_PAGES_BRANCH}."
_info "      If not already configured, set GitHub Pages to serve from that branch (Settings → Pages)."

# End of script
