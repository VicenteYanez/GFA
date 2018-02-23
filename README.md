# GNSS Field Analysis (GFA)

My name is Vicente Yañez and I'm a geologist interested in computer science. GFA is a compilation of python's functions that I made for my internship in [Centro Sismológico Nacional](http://www.sismologia.cl/) (CSN) and rutines I wrote during my undergraduate thesis.

GFA have the functions to calculate vectors from a gnss time series data with the trajectory model of [Bevis & Brown, 2014](https://link.springer.com/article/10.1007/s00190-013-0685-5) and with this, obtain a velocity field that represent the cortical deformation. Also, for the calculation of the velocity field, GFA has the distance-weighed function from [Cardozo & Allmendinger](http://www.sciencedirect.com/science/article/pii/S0098300408002410).After, from the velocity field this program include the tools to calculate the velocity tensor, and from this, the vorticity and stretching tensors (from the equations of [Davis & Titus, 2011](http://www.joshuadavis.us/teaching/2013fcomps/davistitus2011.pdf)).

I am sorry if you read any misspells o grammatical error. English is my second lenguage and my skills are far from perfect.

## Andes 3DB
For the chilean investigation proyect "Active Tectonics and Volcanism at the Southern", a web based GUI was made. The GUI was developed using Python Flask. It provides

![alt text](https://github.com/VicenteYanez/GFA/blob/develop/static/images/homepage.png?raw=true)

![alt text](https://github.com/VicenteYanez/GFA/blob/develop/static/images/homepage.png?raw=true)

## Cientific development using GFA

## Collaborators
Andrés Tassara. Lead investigator of the chilean Fondecyt proyect 1151175 Active Tectonics and Volcanism at the Southern Andes.

Francisco García. Francisco is a geologist doing his doctoral thesis in Universidad de Concepción about the connection of sismicity and volcanic eruptions in the southern Andes. You can contact him to his email franciscogarcia@udec.cl

You can contact me through mail vicenteyanez@protonmail.com

## Install
```
sudo pip3 install virtualenv
mkdir my_project
cd my_project
virtualenv my_project_env
. my_project/bin/activate
sudo pip3 install gfa
```
Through pip install the required libraries automatically, some times I had some problems installing the cartopy library. If that is your case I recomend installing the libraries one by one using pip before install GFA.
### Requirements
* numpy
* scipy
* matplotlib
* cartopy
* click
* flask
* pandas
### For use it on a PC or laptop

### For use it on a server
All the dependences necessary for run GFA are installable via pip. Also, the use of virtualev is recommended.



## Roadmap
