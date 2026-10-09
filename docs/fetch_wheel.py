"""Download the wheel that the ``docs`` workflow built for this commit.

Run by ``.readthedocs.yml`` as ``python docs/fetch_wheel.py <directory>``. It
needs ``GH_ACTIONS_TOKEN`` and the commit Read the Docs checked out
(``READTHEDOCS_GIT_COMMIT_HASH``). It waits for the workflow while it runs, and
exits with status 1 when there is no wheel to use, which stops the docs build.
"""

from __future__ import annotations

import io
import json
import os
import sys
import time
import urllib.error
import urllib.request
import zipfile

API = "https://api.github.com/repos/bhardwaj-lab/sincei"
WORKFLOW = "docs.yml"
ARTIFACT = "docs-wheel"
POLL = 20
TIMEOUT = 10 * 60


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args: object, **kwargs: object) -> None:
        return None


def _request(url: str, token: str) -> urllib.request.Request:
    return urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "sincei-docs",
        },
    )


def _json(url: str, token: str) -> dict:
    with urllib.request.urlopen(_request(url, token), timeout=60) as response:
        return json.load(response)


def _download(url: str, token: str) -> bytes:
    opener = urllib.request.build_opener(_NoRedirect)
    try:
        with opener.open(_request(url, token), timeout=60) as response:
            return response.read()
    except urllib.error.HTTPError as redirect:
        if redirect.code not in (301, 302, 303, 307, 308):
            raise
        location = redirect.headers["Location"]
    with urllib.request.urlopen(location, timeout=300) as response:
        return response.read()


def main(directory: str) -> int:
    token = os.environ.get("GH_ACTIONS_TOKEN")
    commit = os.environ.get("READTHEDOCS_GIT_COMMIT_HASH")
    if not token or not commit:
        print("fetch_wheel: GH_ACTIONS_TOKEN or the commit is not set", file=sys.stderr)
        return 1

    deadline = time.monotonic() + TIMEOUT
    while True:
        runs = _json(
            f"{API}/actions/workflows/{WORKFLOW}/runs?head_sha={commit}&per_page=10",
            token,
        )["workflow_runs"]
        for run in runs:
            artifacts = _json(f"{API}/actions/runs/{run['id']}/artifacts", token)
            for artifact in artifacts["artifacts"]:
                if artifact["name"] == ARTIFACT and not artifact["expired"]:
                    archive = _download(artifact["archive_download_url"], token)
                    with zipfile.ZipFile(io.BytesIO(archive)) as wheels:
                        wheels.extractall(directory)
                        names = ", ".join(wheels.namelist())
                        print(f"fetch_wheel: {names} from run {run['id']}")
                    return 0

        finished = runs and all(run["status"] == "completed" for run in runs)
        if finished or time.monotonic() > deadline:
            break
        print(f"fetch_wheel: waiting for the docs workflow of {commit}")
        time.sleep(POLL)

    print(
        f"fetch_wheel: no {ARTIFACT} artifact for {commit}; did the docs workflow "
        "run for this commit?",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
