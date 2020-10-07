# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey and OneToOneField has `on_delete` set to the desired behavior
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.

from django.db import models
from dal import autocomplete


# Create your models here.

class Experiment(models.Model):
    exid = models.BigAutoField(db_column='exId', primary_key=True)  # Field name made lowercase.
    exidname = models.CharField(db_column='exIdName', max_length=30)  # Field name made lowercase.
    patientid = models.ForeignKey('Patient', models.DO_NOTHING, db_column='patientId')  # Field name made lowercase.
    type = models.CharField(max_length=18)

    def __str__(self):
        return self.exidname
    
    class Meta:
        managed = False
        db_table = 'experiment'


class Geneexpression(models.Model):
    geneid = models.BigAutoField(db_column='geneId', primary_key=True)  # Field name made lowercase.
    genename = models.CharField(db_column='geneName', max_length=30)  # Field name made lowercase.
    geneexpressionvalues = models.FloatField(db_column='geneExpressionValues')  # Field name made lowercase.
    exid = models.ForeignKey(Experiment, models.DO_NOTHING, db_column='exId')  # Field name made lowercase.

    #search_gene = models.ManyToManyField('catalog.genename')
    #widgets = {'search_gene': autocomplete.ModelSelect2Multiple(url='geneautocomplete')}

    def __str__(self):
        return self.genename
    
    class Meta:
        managed = False
        db_table = 'geneExpression'


class Geneset(models.Model):
    genesetid = models.BigAutoField(db_column='genesetId', primary_key=True)  # Field name made lowercase.
    genesetname = models.CharField(db_column='genesetName', max_length=50)  # Field name made lowercase.
    geneid = models.ForeignKey(Geneexpression, models.DO_NOTHING, db_column='geneId')  # Field name made lowercase.

    def __str__(self):
        return self.genesetname
    
    class Meta:
        managed = False
        db_table = 'geneset'


class Patient(models.Model):
    patientid = models.BigAutoField(db_column='patientId', primary_key=True)  # Field name made lowercase.
    paidname = models.CharField(db_column='paIdName', max_length=30)  # Field name made lowercase.
    artype = models.CharField(db_column='arType', max_length=10)  # Field name made lowercase.


    def __str__(self):
        return self.paidname
    
    class Meta:
        managed = False
        db_table = 'patient'

