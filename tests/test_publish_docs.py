"""Tests for ``publish_docs_gh-pages.sh``.

The script publishes a built site to a Pages branch as a single parentless commit,
force-pushed. Everything it does is exercised here against a bare repository used as
the remote, so the tests need no network and no GitHub credentials.

"""

import pathlib
import shutil
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent.parent
SCRIPT = REPO / "publish_docs_gh-pages.sh"
HOMEPAGE_GITIGNORE = REPO / "test_example" / "homepage" / ".gitignore"

# Files written into the site directory by the `site_repo` fixture. The leading
# underscore and the nested directory are both cases the published branch must keep.
SITE_FILES = {
    "index.html": "<html>index</html>",
    "htmls/REGN10933_mut_effect.html": "<html>escape</html>",
    "notebooks/sub/nested.html": "<html>nested</html>",
    "_leading_underscore.html": "<html>underscore</html>",
}

# A VitePress homepage as the pipeline places it in `results/publish_docs`: the source
# tree with `index.md` and no `index.html`, which its own workflow builds on the Pages
# branch. The dot-prefixed directories must be published.
VITEPRESS_FILES = {
    "index.md": "---\nlayout: home\n---\n",
    "cell_entry.md": "# Cell entry\n",
    "package.json": "{}\n",
    ".vitepress/config.mjs": "export default {};\n",
    ".github/workflows/deploy-gh-pages.yaml": "name: Deploy\n",
    "public/appendix.html": "<html>appendix</html>",
    "public/htmls/_leading_underscore.html": "<html>underscore</html>",
}

# Files that `npm install` and `vitepress build` leave in a homepage, all excluded by
# the homepage's own `.gitignore`.
VITEPRESS_BUILD_FILES = {
    "node_modules/vitepress/package.json": "{}\n",
    "package-lock.json": "{}\n",
    ".vitepress/cache/deps/vue.js": "// cache\n",
    ".vitepress/dist/index.html": "<html>built</html>",
}


def git(*args, cwd):
    """Run a `git` command, returning its stripped stdout."""
    return subprocess.run(
        ["git", *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def run_script(cwd, env_extra=None, script=SCRIPT):
    """Run the publish script, returning the `CompletedProcess` without checking it."""
    env = {
        "PATH": "/usr/bin:/bin:/usr/local/bin",
        # Isolate from the invoking user's git configuration, so the test does not
        # depend on their `core.excludesFile`, `init.defaultBranch`, or identity.
        "HOME": str(cwd),
        "GIT_CONFIG_NOSYSTEM": "1",
    }
    env.update(env_extra or {})
    return subprocess.run(
        ["bash", str(script)],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False,  # the tests assert on the return code themselves
    )


@pytest.fixture
def site_repo(tmp_path):
    """A repo with a built site in `results/publish_docs` and a bare repo as remote."""
    remote = tmp_path / "remote.git"
    proj = tmp_path / "proj"
    git("init", "-q", "--bare", str(remote), cwd=tmp_path)
    git("init", "-q", str(proj), cwd=tmp_path)
    git("config", "user.name", "Test User", cwd=proj)
    git("config", "user.email", "test@example.com", cwd=proj)
    git("remote", "add", "origin", str(remote), cwd=proj)

    # The core of the lab default `.gitignore`, which ignores the whole of `results/`
    # and every `_`-prefixed and dot-prefixed file. The published site must contain
    # those files regardless, so this is what would catch the site being staged with
    # the repository's ignore rules in force.
    (proj / ".gitignore").write_text(".*\n_*\n*.out\n!.git*\nresults/*\n")
    git("add", ".gitignore", cwd=proj)
    git("commit", "-q", "-m", "init", cwd=proj)
    git("push", "-q", "origin", "HEAD:refs/heads/main", cwd=proj)
    # `git ls-remote --symref` only reports a default branch if the remote's HEAD
    # resolves, which for a repo created by `git init --bare` it does not until set.
    git("symbolic-ref", "HEAD", "refs/heads/main", cwd=remote)

    write_site(proj, SITE_FILES)
    return proj, remote


def site_dir(proj):
    """The project's site directory."""
    return proj / "results" / "publish_docs"


def write_site(proj, files):
    """Write `files` (a path -> contents mapping) into the project's site directory."""
    for name, contents in files.items():
        path = site_dir(proj) / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)


def exclude_everything_in_callers_git_config(proj):
    """Make the caller's `info/exclude` and `core.excludesFile` each ignore all files."""
    (proj / ".git" / "info" / "exclude").write_text("*\n")
    excludes_file = proj.parent / "global_excludes"
    excludes_file.write_text("*\n")
    git("config", "core.excludesFile", str(excludes_file), cwd=proj)


def published_files(remote, branch="gh-pages"):
    """The set of file paths on `branch` of the remote."""
    listing = git("ls-tree", "-r", "--name-only", branch, cwd=remote)
    return set(listing.splitlines())


def test_publishes_exact_snapshot(site_repo):
    """The branch holds exactly the site's files, plus `.nojekyll`."""
    proj, remote = site_repo
    result = run_script(proj)
    assert result.returncode == 0, result.stderr
    assert published_files(remote) == set(SITE_FILES) | {".nojekyll"}
    for name, contents in SITE_FILES.items():
        assert git("cat-file", "-p", f"gh-pages:{name}", cwd=remote) == contents


def test_snapshot_commit_has_no_history(site_repo):
    """Publishing twice leaves one parentless commit, not a growing history."""
    proj, remote = site_repo
    assert run_script(proj).returncode == 0
    assert run_script(proj).returncode == 0
    assert git("rev-list", "--count", "gh-pages", cwd=remote) == "1"
    assert git("log", "-1", "--format=%P", "gh-pages", cwd=remote) == ""


def test_republishing_drops_deleted_files(site_repo):
    """A file removed from the site directory disappears from the branch."""
    proj, remote = site_repo
    assert run_script(proj).returncode == 0
    (site_dir(proj) / "_leading_underscore.html").unlink()
    write_site(proj, {"index.html": "<html>index v2</html>"})

    assert run_script(proj).returncode == 0
    assert "_leading_underscore.html" not in published_files(remote)
    assert (
        git("cat-file", "-p", "gh-pages:index.html", cwd=remote)
        == "<html>index v2</html>"
    )


def test_publishes_vitepress_source_honoring_its_gitignore(site_repo):
    """A VitePress source tree publishes with its `.gitignore`'d files left out.

    The caller's own excludes ignore everything here, to show they play no part: only
    the site's `.gitignore` decides what is left out.

    """
    proj, remote = site_repo
    shutil.rmtree(site_dir(proj))
    write_site(proj, VITEPRESS_FILES)
    write_site(proj, VITEPRESS_BUILD_FILES)
    shutil.copy(HOMEPAGE_GITIGNORE, site_dir(proj) / ".gitignore")
    exclude_everything_in_callers_git_config(proj)

    result = run_script(proj)
    assert result.returncode == 0, result.stderr
    assert published_files(remote) == set(VITEPRESS_FILES) | {
        ".gitignore",
        ".nojekyll",
    }


def test_leaves_no_trace_in_the_calling_repo(site_repo):
    """No local branch, staging ref, worktree, or scratch survives the publish."""
    proj, _remote = site_repo
    tmpdir = proj.parent / "tmpdir"
    tmpdir.mkdir()

    assert run_script(proj, {"TMPDIR": str(tmpdir)}).returncode == 0

    assert git("branch", "--list", "gh-pages", cwd=proj) == ""
    assert git("for-each-ref", "refs/publish-pages", cwd=proj) == ""
    assert git("worktree", "list", "--porcelain", cwd=proj).count("worktree ") == 1
    assert list((proj / ".git").glob("publish-pages.*")) == []
    # The script must not use TMPDIR at all: a site of any size needs no scratch space.
    assert list(tmpdir.iterdir()) == []


def test_does_not_disturb_another_runs_refs(site_repo):
    """A publish leaves refs under `refs/publish-pages/` alone.

    The snapshot commit is pushed by object name, so nothing is written under that
    namespace and nothing is cleaned up from it. Were a shared ref introduced there,
    together with a cleanup that deletes it unconditionally, two overlapping publishes
    could delete each other's in-flight commit.

    """
    proj, _remote = site_repo
    other = git("rev-parse", "HEAD", cwd=proj)
    git("update-ref", "refs/publish-pages/other-run", other, cwd=proj)

    assert run_script(proj).returncode == 0

    assert git("rev-parse", "refs/publish-pages/other-run", cwd=proj) == other


def test_callers_excludes_cannot_drop_files(site_repo):
    """`$GIT_DIR/info/exclude` and `core.excludesFile` cannot drop published files."""
    proj, remote = site_repo
    exclude_everything_in_callers_git_config(proj)

    assert run_script(proj).returncode == 0
    assert published_files(remote) == set(SITE_FILES) | {".nojekyll"}


def test_refuses_to_publish_a_snapshot_with_no_index(site_repo, tmp_path):
    """A snapshot missing `index.html` and `index.md` is refused, not force-pushed.

    The force-push would otherwise replace a live site with nothing. Staging is
    sabotaged here the only way that can produce an empty snapshot -- making it read
    the caller's excludes, which ignore everything -- to prove the guard is reachable
    and that the published branch survives.

    """
    proj, remote = site_repo
    assert run_script(proj).returncode == 0
    published_before = published_files(remote)

    script_text = SCRIPT.read_text()
    sabotaged_text = script_text.replace(
        "--exclude-per-directory=.gitignore", "--exclude-standard"
    )
    assert sabotaged_text != script_text
    sabotaged = tmp_path / "sabotaged.sh"
    sabotaged.write_text(sabotaged_text)
    exclude_everything_in_callers_git_config(proj)

    result = run_script(proj, script=sabotaged)
    assert result.returncode != 0
    assert "no top-level index.html or index.md" in result.stderr
    assert published_files(remote) == published_before


def test_refuses_the_remotes_default_branch(site_repo):
    """Publishing over the remote's default branch is refused."""
    proj, _ = site_repo
    result = run_script(proj, {"PUBLISH_DOCS_GH_PAGES_BRANCH": "main"})
    assert result.returncode != 0
    assert "default branch 'main'" in result.stderr


@pytest.mark.parametrize("branch", ["master", "develop"])
def test_refuses_reserved_branch_names(site_repo, branch):
    """Publishing to a branch that is conventionally a source branch is refused."""
    proj, _ = site_repo
    result = run_script(proj, {"PUBLISH_DOCS_GH_PAGES_BRANCH": branch})
    assert result.returncode != 0
    assert f"Refusing to publish to '{branch}'" in result.stderr


def test_requires_an_index_in_the_site_dir(site_repo):
    """A site directory with no top-level `index.html` or `index.md` is rejected."""
    proj, _ = site_repo
    (site_dir(proj) / "index.html").unlink()
    result = run_script(proj)
    assert result.returncode != 0
    assert "Missing index.html or index.md" in result.stderr


def test_requires_the_site_dir_to_exist(site_repo):
    """A missing site directory is rejected up front."""
    proj, _ = site_repo
    result = run_script(proj, {"PUBLISH_DOCS_GH_PAGES_SITE_DIR": "results/nosuch"})
    assert result.returncode != 0
    assert "Directory not found" in result.stderr


def test_requires_git_identity(site_repo):
    """Without `user.name`/`user.email` the script stops before touching the remote."""
    proj, _ = site_repo
    git("config", "--unset", "user.email", cwd=proj)
    result = run_script(proj)
    assert result.returncode != 0
    assert "user.email is not configured" in result.stderr
