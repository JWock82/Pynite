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

## Testing

Run tests using `pytest`.

```bash
pytest
```
