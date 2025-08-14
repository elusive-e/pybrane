**Welcome to PyBRANE: Python Based Reasearch/ANalysis Engine (For Computational Chemistry)!**
NEW UPDATE COMING SOON!!!!!! PRE BETA RELEASE EXPECTED AROUND SEPTEMBER!!!

LICENSE: Apache 2. Read more here: https://www.apache.org/licenses/LICENSE-2.0.html

Thanks for checking this out!

**Overview:**
  - GUI with rich functionality*
  - 🌟 Visual outputs of molecules with mouse support*
  - 2D and 3D views of molecules
  - Ability to build bond and atoms with menus
  - 🌟Functional section that allows for the building of a custom biological membrane with backend scripts
  - Error logs and error handling help the program avoid crashing
  - Web file retrival*
  - 🌟Shell section for running commands in the program
  - Simulation tab powered by OpenMM*
  - 🌟Analysis Engine with automatic graphs detailing the results of the experiment*
  - Chat bot to help with inquires*
items marked with * are currently REALLY buggy (pre-beta)
items marked with 🌟 are highlight features

**How this program was created:**
I'm elusive E and im a highschooler from the eastern seaboard. When i went to the science fair last year, i was very intrigued by the molecular simulation projects and how they got such high marks. I was also really confused when one guy next to me was told he needed to do more of the project himself and he had a computational biology project with GROMACS...

ANYWAYS! I tried to set up GROMACS on my machine, but only god knows how it got there, honestly. I looked into a bunch of different tools such as YASARA, VMD, GROMACS, CHARMGUI, etc, but they were either hard to use or they were 2000 dollars. I'm lucky i have a local community college email and got to have access to some of the programs that gave out free uni licenses, because without those i would be double screwed.

The point is, seeing all that made me kinda sad about how something so cool is so freaking hard. So at first this all started as a script to make a more realistic e. coli membrane with embeeded proteins, but then i just kept going. My setup for the project was Thonny because my VS kinda broke. 


AI Disclosure: Some sections of code have been generated with the help of AI. I have tried my best to use AI as a way to have examples for the libraries I need, but some sections are heavy refrenced from AI. I had no coding friends or people who knew about this, so AI helped me get through tough spots.
SECOND AI Disclosure: I AM trying to make a chatbot to help with answering questions relating to the program, however it is quite broken right now because I messed up some dependencies a while back. 

**Known Issues and Future Improvements:**
Currently, the program will crash if theres an error I havet accounted for. I've tried my best to give warnings to the user, but sometimes life happens. This is my first release so I'm not really sure how many bugs are acceptable, but lets say the exterminators would have a field day.

As for improvements, the 3d OpenGL viewer with the molecules is very fussy with large molecules. So I need to go back and optimize the matrix handling, as well as see if there is a more efficient way to export moleulces to the viewer.
Some other improvements:
  - Finish adding in all the little stuff in the top menu bar
  - Add in a window for file conversions for simulations
  - Add in a window for search customiziation
  - Allow for more customization with the graphs
  - Finish the animation feature on the analysis viewer
  - Rename variables so that they aren't so confusing + more consistent

**Michael: The AI**
This is going to be breif, but MICHAEL (Machine Intelligent CHAtbot Engine (for) Learning) is the inhouse AI for answering questions abotu how to use the program and quesitons on molecular dynamics/modeling. He is currently really terrible at doing that. I'm sorry.

**How to Install/Use:**
Im not going to assume someone reading this knows about virtual enviorments, because I know I didn't before this project.

Conda/miniconda works great for me, and allows you to download a lot of fun stuff. 

Link to the official documentation: https://docs.anaconda.com/miniconda/

Python version I wrote in: 3.11.5 (used to be higher, then i ducked up my virtual enviorment and couldn't fix it)

Once you get miniconda, install the requirments in requirments.txt, and you should be good to go. Run the main.py script, I've gone through several migranes with passing functions back and forth, so everything should get set up automatically.

If you run the ui set up file by accident, you will have no backend features, so make sure you're running main

I also HIGHLY reccomend you keep the program folder somewhere where you knwo how to find it, as all files are currently saved into that folder.

**How to Contribute:**
If you like what I'm doing, and you want to help me make the project better, I'm eternally grateful. Just using and downloading it is already intensely satifying, but if you want to use the code or help make it run smoother, do that. 

**More about me and my story:**
feel free to reach me at imelusivee@gmail.com for any specific requests. if you have a similar tool, i would love to hear more about it! please hit me up

As mentioned, I'm a highschooler and i made this project over the summer. Now that school is starting i still intend to work on it, however the drive has somewhat left me after i've created mostly what i wanted to. What was my motivation? Well I want to help make MD and MM a bit more accesable with this tool to help people get their toes wet, so they have enough confidence to get into the really complicated and powerful tools. I hope people find this helpful, in some way.

I did this project solo, and since i didn't really have anyone to help me i turned to ai for help on topics i couldn't understand. I'm not proud of this fact, but it is a fact and its important for me to be transparent. I've been commited to understanding any code i've been putting in from outside sources, not limited to AI. I can explain how it works, and i know what each section does. I'm not close to finishing highschool so some math and other topics for understanding some code is still out of reach for me, but I've tried my best.

This also applies to the field of computational chem, which i decided to make a tool for. I don't really know if theres a dire need for it, and maybe someone will look at this and be like "bro are you 
stupid we have that already people dont struggle with this" but oh well, I can't take back those hours now. I've finished geometry and bio, taking ap stat and chem this year, so i don't know everything about molecular simualtions and how they work. But i'm trying.

Why did I choose to publish now? Well, now that i consider the project in beta pretty much I want to see if people use it, or if i wasted my time lmao (sobs). I also think now that the initial setup is done, I'm ready to open myself upto critique because that's important if i want to move forward. Before, i had more key features i needed to add, but now its a lot of polishing (which i suck at).


All in all, thanks for reading if you made it this far, here are some screenshots of ICONIC moments of the AI chatbot, Michael:
*pictures will appear soon*

