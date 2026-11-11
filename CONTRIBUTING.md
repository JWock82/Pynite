If you would like to contribute to Pynite, please review the following guidelines.

## Coding Guidelines:
1. Over-commenting is better than under commenting. Try to provide a comment for nearly every line of code.
2. Keep it organized. Large methods (functions/subroutines) can often be broken into smaller pieces. Classes and inheritance can also reduce unecessary repetition. Place things where they make sense, and provide method and variable names that are descriptive.
3. Keep it clean. Try to follow the PEP8 style guide as much as possible. There are some conventions in there that don't work well for finite element analysis, where we use capitalization to distinguish between local and global element terms, but otherwise please stick to the guide. I also dislike the limitation on line widths of 79 characters. I end up breaking a lot of lines if I do that. I'm using 100 characters. In some instances I use more (such as for matrix rows that really need to be displayed on one line to make sense). It's recognized that Pynite already has room for improvement in this regard, but going forward PEP8 is the goal.

## Keys to Getting Your Contributions Integrated Quickly:
1. First and foremost, **KEEP IT SIMPLE**. Every complex problem is nothing more than a series of small simple problems building on each other. Try to break your code into small self-explanatory pieces. Consider issuing larger features as a series of smaller features that add incremental functionality that builds on prior functionality, rather than overhauling the program radically all at once. The addition of CI to Pynite has helped check new code for unexpected issues, but I still like to review this stuff before accepting it. Pull-requests get approved much faster if they are easy for me to follow.
2. Provide features that most users would find useful.
3. Test and debug your code. Please make sure it's working properly before submitting a pull request.
4. Provide an example problem `.py` file with your code.
5. Provide all pieces necessary in your pull request.
6. Follow the coding guidelines above. This will make review go much faster.
7. Pull requests for any active repository `Projects` will usually be given first priorty.

## AI-Assisted Contributions

Recently, Pynite has been receiving more contributions than it has in the past. AI is assisting with many if not most of these contributions. AI tools are helpful and may be used to help develop contributions, but it should be used with caution. Pynite is structural engineering software, and AI generated code should be carefully reviewed by humans, line by line, with expertise.

For AI-assisted pull requests:

1. **Disclose AI assistance.** In the pull request description, identify the tools used and briefly explain which parts of the contribution they helped produce or review. Routine editor autocomplete does not need to be disclosed.
2. **Understand the engineering behind the contribution.** AI is not a substitute for human engineering knowledge or judgment. Where a change involves engineering, understand its underlying principles, assumptions, limitations, and intended use well enough to explain and evaluate them yourself. Do not rely on AI-generated explanations to defend engineering work you do not understand.
3. **Review and understand the whole change.** Read every changed line, including tests and documentation. Be prepared to explain the implementation, assumptions, and design choices, and to respond to review feedback. Do not submit a pull request you cannot maintain or explain. If you are using Visual Stuio Code, consider setting your session to a "Local" harness session to allow for line-by-line review before accepting sweeping changes.
4. **Avoid bloat.** AI tools often add unnecessary code, abstractions, dependencies, generated files, or unrelated changes. Remove them before submitting. Keep the contribution as small and simple as the problem allows.
5. **Keep the code understandable.** AI-generated code can obscure its intent behind complicated logic, indirection, or vague names. Simplify it, use descriptive names, and make the implementation easy for a maintainer to follow. Do not submit code you cannot explain.
6. **Test the actual change.** Test the change yourself. Don't rely on AI written tests unless you have reviewed them and find them to be meaningful and adequate.
7. **Comment thoroughly.** AI tools often omit comments. Follow the coding guidelines above: aim for a comment on nearly every line, and explain the reasoning, equations, and assumptions behind the implementation. Review AI generated comments for accuracy and usefulness.
8. **Check provenance and confidentiality.** Submit only material you have the right to contribute under this project's license. Do not put confidential or third-party material into an AI tool without permission, and verify any references or citations before including them.

Maintainers may ask for clarification or decline contributions that are bloated, difficult to understand, inadequately commented, or not sufficiently reviewed or verified, regardless of how they were produced.

## Development Setup

Pynite uses [uv](https://docs.astral.sh/uv/) to manage development environments and locked
dependencies. Install uv, then create the environment and install the project with its test tools
and optional runtime dependencies:

```bash
uv sync --locked --extra all --all-groups
```

When project dependencies change, update `pyproject.toml` and regenerate the committed lockfile:

```bash
uv lock
```

## Testing

Run tests inside the managed environment using `uv run`:

```bash
uv run --locked pytest
```

## Linting

Linting is handled by [ruff](https://docs.astral.sh/ruff/), with hooks managed by [prek](https://github.com/j178/prek). The rules live in `ruff.toml`.

Install prek outside the project environment, so that git GUI clients can run the hooks without activating a virtualenv:

```bash
uv tool install prek
```

`pipx install prek` works too. Then install the git hook:

```bash
prek install
```

The configured checks now run against your staged files on every commit, fixing what they can. To run them across the whole repository:

```bash
prek run --all-files
```

The same checks run in CI on every pull request.

prek reads `.pre-commit-config.yaml`, which is the standard pre-commit config format, so stock `pre-commit` works against the same config if you already use it. If you want ruff in your editor, install it separately and match the version pinned in `.pre-commit-config.yaml` so it reports the same results as the hook.
