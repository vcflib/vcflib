This is a dead simple progress bar implementation for C. 
It has no extra dependencies (like libtermcap / ncurses) 
and should be easily dropped into any project. As an added
bonus, it'll report time ETAs for long running jobs.


The main interface is `print_progress`, which takes in a percentage
and a timestamp. Call it repeatedly in a loop to make the progress
bar grow across the screen.

![](example.gif)
