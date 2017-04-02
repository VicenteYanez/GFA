# GNSS Field Analysis (GFA)

My name is Vicente Yañez and I'm a geologist interested in computer science. GFA is a compilation of python's functions that I wrote for my internship in [Centro Sismológico Nacional](http://www.sismologia.cl/) (CSN) and during my undergraduate thesis.

GFA have the functions to calculate vectors from a gnss time series data with the trajectory model of [Bevis & Brown, 2014](https://link.springer.com/article/10.1007/s00190-013-0685-5) and with this, obtain a velocity field that represent the cortical deformation. After, from the velocity field this program include the tools to calculate the velocity tensor, and from this, the vorticity and stretching tensors (from the equations of [Davis & Titus, 2011](http://www.joshuadavis.us/teaching/2013fcomps/davistitus2011.pdf)).

To start using GFA, you should see the wiki page [How to start](https://github.com/VicenteYanez/GFA/wiki/How-to-start)

I'm sorry if you read any misspells o grammatical error. English is my second lenguage and my skills in it is far from perfect.

## Collaborators
Francisco García. Francisco is a geologist who is doing his doctoral thesis in Universidad de Concepción about the connection of sismicity and volcanic eruptions in the southern Andes. You can contact him to his email franciscogarcia@udec.cl

Also you can contact me to my mail vicenteyanez@protonmail.com

# Install

All the dependences necessary for run GFA are already contained in gfa via virtualenv or optionally installable with pip.

```
git clone http://github.com/VicenteYanez/GFA.github
cd GFA
. gfa_env/bin/activate
```
or if you want to install the libraries on your computer
```
sudo pip3 install -r requirements.txt
```

If everything is alright, you can continue in the [How to start page](https://github.com/VicenteYanez/GFA/wiki/How-to-start)
