from django.db import models

# Create your models here.
class Strain(models.Model):
    name = models.CharField(max_length = 200)
    taxonid = models.CharField(max_length = 10)
    pmid = models.CharField(max_length = 20)
    createdDate = models.DateTimeField('created date')
    modifiedDate = models.DateTimeField('modified date')

    def __unicode__(self):
        return self.name
    # Strain ##


class Contig(models.Model):
    strain = models.ForeignKey(Strain)
    name = models.CharField(max_length = 500)
    seq = models.TextField()
    createdDate = models.DateTimeField('created date')
    modifiedDate = models.DateTimeField('modified date')

    def __unicode__(self):
        return self.name

    def getSize(self):
        return len(self.seq)
    getSize.short_description = 'Contig Length'
    ## Contig ##

class Reference(models.Model):
    sgdid = models.CharField(max_length = 15)
    feature_type = models.CharField(max_length = 100)
    qualifier = models.CharField(max_length = 100, null = True)
    feature_name = models.CharField(max_length = 200, null = True)
    standard_name = models.CharField(max_length = 200, null = True)
    parent_name = models.CharField(max_length = 200, null = True)
    aliases = models.CharField(max_length = 1000, null = True)
    secondary_sgdid = models.CharField(max_length = 1000, null = True)
    description = models.CharField(max_length = 1000, null = True)
    createdDate = models.DateTimeField('created date')
    modifiedDate = models.DateTimeField('modified date')

class Feature(models.Model):
    reference = models.ForeignKey(Reference, null = True)
    contig = models.ForeignKey(Contig)
    feature_id = models.CharField(max_length = 100, null = True)
    parent = models.ForeignKey('self', null = True)
    feature_type = models.CharField(max_length = 100)
    start_coord = models.IntegerField()
    stop_coord = models.IntegerField()
    strand = models.CharField(max_length = 1)
    createdDate = models.DateTimeField('created date')
    modifiedDate = models.DateTimeField('modified date')
    ## Feature ##

class CDS(models.Model):
    feature = models.ForeignKey(Feature)
    name = models.CharField(max_length = 500)
    seq = models.TextField()
    createdDate = models.DateTimeField('created date')
    modifiedDate = models.DateTimeField('modified date')
    ## CDS ## 

class Protein(models.Model):
    cds = models.ForeignKey(CDS)
    name = models.CharField(max_length = 500)
    seq = models.TextField()
    createdDate = models.DateTimeField('created date')
    modifiedDate = models.DateTimeField('modified date')
    ## Protein ##
