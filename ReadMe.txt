Files needed to modify
catalog folder:
static/css/
styles.css: decorate website

templates folder:
base_generic.html: basic framework
index.html: home page
plot.html: plotting page

models.py: the codes are automatically produced when connected with MySQL database successfully. 

(There are two urls.py files. This one is in catalog folder. The other one is in locallibrary folder)
urls.py: saving the function names and setting url address for corresponding function. 

views.py: saving all of the back-end functions.

locallibrary folder:
settings.py: two parts need to add codes. I put #Adding codes in front of the codes. One is  the INSTALLED_APPS, the other part is in the end of the file.

urls.py: Adding some codes when setting up the website framework.
