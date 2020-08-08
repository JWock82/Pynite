If you would like to contribute to PyNite, please review the following guidelines.

## Coding Guidelines:
1. Keep it simple. Try to keep your code as simple and self-explanatory as possible. Over-commenting is better than under commenting. Try to provide a comment for nearly every line of code.
2. Keep it organized. Spaghetti code is hard to follow. Place things where they make sense, and provide method and variable names that are descriptive.
3. Keep it clean. Try to follow the PEP8 style guide as much as possible. There are some conventions in there that don't work well for finite element analysis, such as capitalization rules for method names, but otherwise please stick to the guide. I also dislike the limitation on line widths of 79 characters. I end up breaking a lot of lines if I do that. I'm using 100 characters. In some instances I use more (such as for matrix rows that really need to be displayed on one line to make sense). It's recognized that PyNite already has room for improvement in this regard, but going forward PEP8 is the goal.

## Keys to Getting Your Contributions Integrated Quickly:
1. Test and debug your code. Please make sure it's working properly before submitting a pull request.
2. Provide an example problem `.py` file with your code.
3. Provide all pieces necessary in your pull request.
4. Keep pull requests limited to one basic feature at a time.
5. Follow the coding guidelines above. This will make review go much faster.
6. Pull requests for any active repository `Projects` will usually be given first priorty.
