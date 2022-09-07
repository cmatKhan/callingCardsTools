# Using tests to make debugging __easy__ and __reproducible__

What I am going to describe has some name -- test first development or some such. 
I don't know what it is called, and I don't care to have a debate about whether 
or not it is "good" or "an effective part of the development cycle". These are 
discussions that happen, apparently, among the "developer community", and to the 
extent that the discussion is interesting, the setting in which they are occuring 
are not academic labs. In other words, if you're tempted to call what I am about 
to describe as "test first development", please set that out of your mind.  

## Tests as debugging tools

I am not suggesting that you write tests that prove that your code is correct, 
nor am I suggesting that you need a certain "test coverage". I am suggesting that 
you use the tools available in your IDE and the language to make writing easier. 
You will be debugging your code, one way or another. If you do so in a test 
environment, then your debugging will be both __easier__ and __reproducible__. 
Why does __reproducible__ matter? It makes it easier to change/update/improve 
your code in the future, and it helps someone else understand what your code 
does without you telling them.  

## Using tests as writing aids

### First, think

Start by thinking about what you want to write. I do this on paper. 

### Write a skeleton

This is an example of a class, but this could also be a set of functions.

### Get some test data

This doesn't have to be perfect. It is something that can be iteratively improved. 
Do be careful about adding/committing large files to git, though -- subset them 
down before adding.

### If you can, in the test, write what you expect to get out of your function

But this can be really hard -- what if you don't know what you are doing yet? 
I am frequently in this position. Best practice would probably be to go back 
to the paper and figure it out, but frequently it is easiest to have an interactive 
environment to play with the data.

### Write a "dummy test", write some code into your function and set a break point

Now you can get into the class/function namespace and start exploring, interactively, 
while you write. You can go back and forth between writing the class/function 
and the test, refining both as you go.  

### But writing a test increases the number of lines that I have to write

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