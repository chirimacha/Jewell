# Bandit
#Instructions for setting environment variable (Mac)
Note: we are currently using "SPATIAL_UNCERTAINTY" as our environment variable 

1. Open Terminal
2. File--> New Command (Shift+Apple+N)
3. In the console that comes up, type:
     launchctl setenv SPATIAL_UNCERTAINTY "....../Bandit" 
 *  "....../Bandit" is the path to the Bandit repository on your computer from GitHub.  
 * NOTE: This path cannot contain any spaces or it will not work. 
4. Run. 
5. You can check to make sure it worked by running a new command: 
launchctl getenv SPATIAL_UNCERTAINTY
 * This should spit out the path to github repository. 
