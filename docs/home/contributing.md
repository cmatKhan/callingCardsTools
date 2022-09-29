
# Contributing

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

## Environment setup

Nothing easier!

Fork and clone the repository, then:

```bash
cd python
make setup
```

> NOTE:
> If it fails for some reason,
> you'll need to install
> [PDM](https://github.com/pdm-project/pdm)
> manually.
> 
> You can install it with:
> 
> ```bash
> python3 -m pip install --user pipx
> pipx install pdm
> ```
> 
> Now you can try running `make setup` again,
> or simply `pdm install`.

You now have the dependencies installed.

You can run the application with `pdm run mkdocstrings-python [ARGS...]`.

Run `make help` to see all the available actions!

## Tasks

This project uses [duty](https://github.com/pawamoy/duty) to run tasks.
A Makefile is also provided. The Makefile will try to run certain tasks
on multiple Python versions. If for some reason you don't want to run the task
on multiple Python versions, you can do one of the following:

1. `export PYTHON_VERSIONS= `: this will run the task
   with only the current Python version
2. run the task directly with `pdm run duty TASK`

The Makefile detects if a virtual environment is activated,
so `make` will work the same with the virtualenv activated or not.

## Development

As usual:

1. create a new branch: `git checkout -b feature-or-bugfix-name`
1. edit the code and/or the documentation

**Before committing:**

1. run `make format` to auto-format the code
1. run `make check` to check everything (fix any warning)
1. run `make test` to run the tests (fix any issue)
1. if you updated the documentation or the project dependencies:
    1. run `make docs-serve`
    1. go to http://localhost:8000 and check that everything looks good
1. follow our [commit message convention](#commit-message-convention)

If you are unsure about how to fix or ignore a warning,
just let the continuous integration fail,
and we will help you during review.

Don't bother updating the changelog, we will take care of this.

## Commit message convention

Commits messages must follow the
[Angular style](https://gist.github.com/stephenparish/9941e89d80e2bc58a153#format-of-the-commit-message):

```
<type>[(scope)]: Subject

[Body]
```

Scope and body are optional. Type can be:

- `build`: About packaging, building wheels, etc.
- `chore`: About packaging or repo/files management.
- `ci`: About Continuous Integration.
- `docs`: About documentation.
- `feat`: New feature.
- `fix`: Bug fix.
- `perf`: About performance.
- `refactor`: Changes which are not features nor bug fixes.
- `style`: A change in code style/format.
- `tests`: About tests.

**Subject (and body) must be valid Markdown.**
If you write a body, please add issues references at the end:

```
Body.

References: #10, #11.
Fixes #15.
```

## Pull requests guidelines

Link to any related issue in the Pull Request message.

During review, we recommend using fixups:

```bash
# SHA is the SHA of the commit you want to fix
git commit --fixup=SHA
```

Once all the changes are approved, you can squash your commits:

```bash
git rebase -i --autosquash master
```

And force-push:

```bash
git push -f
```

If this seems all too complicated, you can push or force-push each new commit,
and we will squash them ourselves if needed, before merging.

## Using tests to make debugging __easy__ and __reproducible__

What I am going to describe has some name -- test first development or some such. 
I don't know what it is called, and I don't care to have a debate about whether 
or not it is "good" or "an effective part of the development cycle". These are 
discussions that happen, apparently, among the "developer community", and to the 
extent that the discussion is interesting, the setting in which they are occuring 
are not academic labs. In other words, if you're tempted to call what I am about 
to describe as "test first development", please set that out of your mind.  

### Tests as debugging tools

I am not suggesting that you write tests that prove that your code is correct, 
nor am I suggesting that you need a certain "test coverage". I am suggesting that 
you use the tools available in your IDE and the language to make writing easier. 
You will be debugging your code, one way or another. If you do so in a test 
environment, then your debugging will be both __easier__ and __reproducible__. 
Why does __reproducible__ matter? It makes it easier to change/update/improve 
your code in the future, and it helps someone else understand what your code 
does without you telling them.  

### Using tests as writing aids

#### First, think

Start by thinking about what you want to write. I do this on paper. 

#### Write a skeleton

This is an example of a class, but this could also be a set of functions.

#### Get some test data

This doesn't have to be perfect. It is something that can be iteratively improved. 
Do be careful about adding/committing large files to git, though -- subset them 
down before adding.

#### If you can, in the test, write what you expect to get out of your function

But this can be really hard -- what if you don't know what you are doing yet? 
I am frequently in this position. Best practice would probably be to go back 
to the paper and figure it out, but frequently it is easiest to have an interactive 
environment to play with the data.

#### Write a "dummy test", write some code into your function and set a break point

Now you can get into the class/function namespace and start exploring, interactively, 
while you write. You can go back and forth between writing the class/function 
and the test, refining both as you go.  

#### But writing a test increases the number of lines that I have to write

That's true, although likely not by much. But, what you have gained is an 
interactive and reproducible debugging and development environment. 
Not only that, anyone else who might want to use your code can do the same. 
It is much easier to have ready made functions (the tests) which bring in some 
data and run your code in a way that allows another person to set a breakpoint 
and take a look around at what your doing. Tests are not specification, and 
they don't take the place of comments and documentation. They don't "prove" 
that your code is right necessarily, and you don't need to spend a ton of time 
coming up with all of the edge cases, etc. The are __tools__ -- they are better 
than print statements for debugging, and they are reproducible in both your 
hands, and the hands of others.

## Cite
This is based on the [mkdocstrings](https://github.com/mkdocstrings/python/blob/master/CONTRIBUTING.md) contribution policy