# GNSS Field Analysis (GFA)

My name is Vicente Yañez and I'm a geologist interested in geophysics and computer science. Right now I am doing a Master in Geophysics in the University of Chile.  The core of GFA is a compilation of python's functions that I made for my internship in [Centro Sismológico Nacional](http://www.sismologia.cl/) (CSN) and rutines I wrote during my undergraduate thesis.You can contact me through mail vicenteyanez@protonmail.com

GFA has the functions to calculate vectors from a gnss time series data with the trajectory model of [Bevis & Brown, 2014](https://link.springer.com/article/10.1007/s00190-013-0685-5) and with this, obtain a velocity field that represent the cortical deformation. Also, for the calculation of the velocity field, GFA has the distance-weighed function from [Cardozo & Allmendinger](http://www.sciencedirect.com/science/article/pii/S0098300408002410).After, from the velocity field this program include the tools to calculate the velocity tensor, and from this, the vorticity and stretching tensors (from the equations of [Davis & Titus, 2011](http://www.joshuadavis.us/teaching/2013fcomps/davistitus2011.pdf)).

The GNSS stations time series Data are not included in GFA.

I am sorry if you read any misspells o grammatical error. English is my second lenguage and my skills are far from perfect.

## Principal commands and first steps
After install GFA, you will have access to the console commands of GFA, you can see them if you write:

```
gfa --help
```
### Set up your config and directories: parameters.ini file
```
gfa configure
```
### GFA automatic solution for all the station in the db
```
gfa buildmodelall
```
### Select your stations
```
gfa select
```
### Plot stations
```
gfa plot
```
### Modify the trajectory model
```
gfa build
```
### Calculate a velocity vector
```
gfa vector
```

### Velocity field analysis

Through there is no command line command to calculate the vorticity or stretching field, you can use the python functions in gfa.field_analysis.field module. Or also you can use the GUI I am making (see next section)
## Andes 3DB
For the chilean investigation proyect "Active Tectonics and Volcanism at the Southern", a web based GUI was made. The GUI was developed using Python Flask. It provides an interface to
1. Select the GNSS station by position
2. Edit and calculate a trajectory model for each station
3. Calculate velocity vectors for differents time intervals.
4. Use this velocity vectors to calculate a velocity field and from this, visualizate the vorticity and streching field.

Screenshot of the selection section. Here you can select the GNSS station you want to use.
![alt text](https://github.com/VicenteYanez/GFA/blob/develop/static/images/homepage1.png?raw=true)

Screenshot of the edition section of the trajectory model parameters and time
![alt text](https://github.com/VicenteYanez/GFA/blob/develop/static/images/homepage2.png?raw=true)

Screenshot of the vorticity visualization on the map
![alt text](https://github.com/VicenteYanez/GFA/blob/develop/static/images/homepage3.png?raw=true)

## Collaborators
Andrés Tassara. Lead investigator of the chilean Fondecyt proyect 1151175 Active Tectonics and Volcanism at the Southern Andes. You can contact him to his mail andrestassara@udec.cl

Francisco García. Francisco is a geologist doing his doctoral thesis in Universidad de Concepción about the connection of sismicity and volcanic eruptions in the southern Andes. You can contact him to his mail franciscogarcia@udec.cl

You can contact me through mail vicenteyanez@protonmail.com

## Install

All the dependences necessary for run GFA are installable via pip. The use of virtualev is recommended.

```
sudo pip3 install virtualenv
mkdir my_project
cd my_project
virtualenv my_project_env
. my_project/bin/activate
sudo pip3 install gfa
```
Through pip install the required libraries automatically, sometimes I had some problems installing the cartopy library. If that is your case I recomend installing the libraries one by one using pip before install GFA.
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



## Roadmap
* Command to calculate vorticity from the command line
* Add fourier spectre analysis to improve the calculation of the automatic solution
* Install and configuration script for its use on the server/local pc
* Real-time edition and visualization of the trajectory model
