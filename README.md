# GNSS Field Analysis (GFA)

My name is Vicente Yañez and I'm a geologist interested in computer science. GFA is a compilation of python's functions that I wrote for my internship in [Centro Sismológico Nacional](http://www.sismologia.cl/) (CSN), during my undergraduate thesis and that I am still developing it.

GFA have the functions to calculate vectors from a gnss time series data with the trajectory model of [Bevis & Brown, 2014](https://link.springer.com/article/10.1007/s00190-013-0685-5) and with this, obtain a velocity field that represent the cortical deformation. Also, for the calculation of the velocity field, GFA has the distance-weighed function from [Cardozo & Allmendinger](http://www.sciencedirect.com/science/article/pii/S0098300408002410).After, from the velocity field this program include the tools to calculate the velocity tensor, and from this, the vorticity and stretching tensors (from the equations of [Davis & Titus, 2011](http://www.joshuadavis.us/teaching/2013fcomps/davistitus2011.pdf)).

I am sorry if you read any misspells o grammatical error. English is my second lenguage and my skills in it is far from perfect.

## Collaborators
Francisco García. Francisco is a geologist who is doing his doctoral thesis in Universidad de Concepción about the connection of sismicity and volcanic eruptions in the southern Andes. You can contact him to his email franciscogarcia@udec.cl

Also you can contact me to my mail vicenteyanez@protonmail.com

# Install

All the dependences necessary for run GFA are installable via pip. Also, the use of virtualev is recommended.

```
sudo pip3 install virtualenv (if you don't have it installed)
mkdir your_project
cd your_project
virtualenv your_project
. your_project/bin/activate
sudo pip3 install gfa
```

If everything is alright, you can continue in the [How to start page](https://github.com/VicenteYanez/GFA/wiki/How-to-start)
