from strains.models import Strain
from strains.models import Contig
from django.contrib import admin

class StrainAdmin(admin.ModelAdmin):
    #fields = ['name', 'taxonid', 'pmid', 'createdDate', 'modifiedDate']
    fieldsets = [ (None,               {'fields' : ['name', 'taxonid', 'pmid']}),
                  ('Date Information', {'fields' : ['createdDate', 'modifiedDate'], 'classes' : ['collapse']}), ]
    list_display = ['name', 'taxonid', 'pmid']

admin.site.register(Strain, StrainAdmin)

class ContigAdmin(admin.ModelAdmin):
    fieldsets = [ (None,               {"fields" : ['strain', 'name', 'seq']}),
                  ('Date Information', {"fields" : ['createdDate', 'modifiedDate'], 'classes' : ['collapse']}), ]
    list_display = ('name', 'getSize')

admin.site.register(Contig, ContigAdmin)

