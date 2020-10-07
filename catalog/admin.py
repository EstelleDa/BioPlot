from django.contrib import admin
from catalog.models import Patient, Experiment, Geneexpression, Geneset

# Register your models here.

class PatAdmin(admin.ModelAdmin):

    list_display = ('patientid', 'paidname', 'artype')
    ordering = ('patientid', )


class ExpAdmin(admin.ModelAdmin):

    list_display = ('exid', 'exidname', 'type')
    ordering = ('exid', )
    search_fields = ('exidname', 'type')

class GeneExpAdmin(admin.ModelAdmin):

    list_display = ('geneid', 'genename', 'geneexpressionvalues')
    ordering = ('geneid', )
    search_fields = ('geneid', 'genename')

class GenesetAdmin(admin.ModelAdmin):

    list_display = ('genesetid', 'genesetname')
    ordering = ('genesetid', )
    search_fields = ('genesetid', 'genesetname')


admin.site.register(Patient, PatAdmin)
admin.site.register(Experiment, ExpAdmin)
admin.site.register(Geneexpression, GeneExpAdmin)
admin.site.register(Geneset, GenesetAdmin)
