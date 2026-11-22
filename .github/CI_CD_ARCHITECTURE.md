# siRNAforge CI/CD architecture

GitHub Actions uses the repository's root `pixi.lock` for Python, native tools,
documentation, and the Linux x86_64 runtime. Locked installs reject manifest
drift. The wheel also supports installation with standard Python tooling.

## Validation workflows

| Workflow | Trigger | Checks |
| --- | --- | --- |
| [CI](workflows/ci.yml) | Pushes and pull requests to `master` or `dev`, or called by the release workflow | Python 3.10–3.12 smoke tests, release host selection, lint, audit, docs, wheel and sdist verification, isolated wheel install and offline container tests |
| [Release candidate](workflows/release-candidate.yml) | Pull requests to `master` and manual dispatch | Version and changelog identity, packaging metadata and the full `make test-release` gate |
| [Development image](workflows/docker-dev.yml) | Pushes to `dev`, matching pull requests to `dev`, and manual dispatch | Builds and tests the locked image. Pushes the tested image on runs outside pull requests. |
| [Documentation](workflows/docs.yml) | Pushes and pull requests to `master`, and manual dispatch | Builds Sphinx from the locked docs environment. The publishing steps apply their branch and event conditions. |

The reusable CI workflow exposes the `CI Status Summary` result to release jobs.
Live public API checks are local opt-ins through `make test-requires-network`.

## Tagged releases

[The release workflow](workflows/release.yml) runs for `v*` tags. It checks tag
and package identity, the dated changelog entry, and ancestry from `master`.
After CI passes, the build job prepares a verified wheel and sdist. The image
build and release tests consume that distribution artifact and check the wheel
digest and source fingerprint.

After the release gate, the workflow creates a draft GitHub release and publishes
the verified versioned artifacts. It checks the published wheel hash and image
identity before making the GitHub release public. The `:latest` image advances
after successful verification of a stable release.

## Local commands

```bash
make dev
make lint
make test-ci
make test-release
```

The full release gate requires Docker. Its container suite runs with networking
disabled. See the [development guide](../docs/developer/development.md) for
locked setup, dependency updates, and rollback.
