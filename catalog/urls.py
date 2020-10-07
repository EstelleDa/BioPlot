from django.conf.urls import url, include
from django.urls import path
from catalog import views
from catalog.views import GeneNameAutocomplete

urlpatterns = [
    path('', views.index, name='index'),
    path('plot/', views.plot, name='plot'),
    path('subplot/', views.subplot, name='subplot'),
    path('removegene/', views.removeGeneFromPlot, name='removeGene'),
    path('addgene/', views.addGeneToPlot, name='addGene'),
    path('removesample/', views.removeSampleFromPlot, name='removeSample'),
    path('addsample/', views.addSampleToPlot, name='addSample'),
]
