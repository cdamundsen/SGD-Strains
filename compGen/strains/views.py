from django.http import HttpResponse, HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.template import Context, loader, RequestContext
from django.shortcuts import render_to_response, get_object_or_404

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalwCommandline


from models import CDS
from models import Contig
from models import Feature
from models import Protein
from models import Reference
from models import Strain

import base64
import bz2
import datetime
import os
import tempfile
import urllib


def index(request):
    strains = Strain.objects.all().order_by('name')
    referenceCount = len(Reference.objects.all())
    
    
    strainList = []
    for strain in strains:
        contigs = Contig.objects.filter(strain = strain)
        features = Feature.objects.filter(contig__strain = strain)
        codingSeqs = CDS.objects.filter(feature__contig__strain = strain)
        proteins = Protein.objects.filter(cds__feature__contig__strain = strain)
        strainList.append([len(contigs), [strain, len(contigs), len(features), len(codingSeqs), len(proteins)]])
    strainList.sort()
    strainList = [x[1] for x in strainList]
    t = loader.get_template('strains/index.html')
    c = Context({'strain_list' : strainList,
                 'reference_count' : referenceCount })
    return HttpResponse(t.render(c))
    ## index ##

def load_references(request):
    return render_to_response('strains/load_references.html',
                              {},
                              context_instance = RequestContext(request))
    ## load_references ##


def save_references(request):
    """
    SGD_features.tab format looks like this:
         0) Primary SGDID (mandatory)
         1) Feature type (mandatory)
         2) Feature qualifier (optional)
         3) Feature name (optional)
         4) Standard gene name (optional)
         5) Alias (optional, multiples separated by |)
         6) Parent feature name (optional)
         7) Secondary SGDID (optional, multiples separated by |)
         8) Chromosome (optional)
         9) Start_coordinate (optional)
        10) Stop_coordinate (optional)
        11) Strand (optional)
        12) Genetic position (optional)
        13) Coordinate version (optional)
        14) Sequence version (optional)
        15) Description (optional)
    """
    featuresFile = request.FILES['features_file']
    features = [x.split('\t') for x in featuresFile.read().splitlines()]
    for line in features:
        reference = Reference()
        reference.sgdid = line[0]
        reference.feature_type = line[1]
        if line[2]:
            reference.qualifier = line[2]
        if line[3]:
            reference.feature_name = line[3]
        if line[4]:
            reference.standard_name = line[4]
        if line[5]:
            reference.aliases = line[5]
        if line[6]:
            reference.parent_name = line[6]
        if line[7]:
            reference.secondary_sgdid = line[7]
        if line[15]:
            reference.description = line[15]

        reference.createdDate = datetime.datetime.now()
        reference.modifiedDate = datetime.datetime.now()

        reference.save()
    return HttpResponseRedirect('/strains/')
    ## save_references ##
    

def detail(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    contigList = Contig.objects.filter(strain = strain_id)
    contigCount = len(contigList)
    featureList = Feature.objects.filter(contig__strain__pk = strain_id)
    featureCount = len(featureList)
    cdsList = CDS.objects.filter(feature__contig__strain__pk = strain_id)
    cdsCount = len(cdsList)
    pepList = Protein.objects.filter(cds__feature__contig__strain__pk = strain_id)
    pepCount = len(pepList)

    t = loader.get_template('strains/detail.html')
    c = Context({ 'the_strain' : theStrain,
                  'contig_count'  : contigCount,
                  'feature_count' : featureCount,
                  'cds_count'     : cdsCount,
                  'pep_count'     : pepCount })

    return HttpResponse(t.render(c))
    ## detail ##

def new_strain(request):
    return render_to_response('strains/new_strain.html',
                              {},
                              context_instance = RequestContext(request))
    ## new_strain ##


def save_new_strain(request):
    strain_name = request.POST['strain_name']
    taxon_id = request.POST['taxon_id']
    pmid = request.POST['pmid']

    s = Strain()
    if strain_name:
        s.name = strain_name
        s.createdDate = datetime.datetime.now()
        s.modifiedDate = datetime.datetime.now()
    if taxon_id:
        s.taxonid = taxon_id
    if pmid:
        s.pmid = pmid
    if strain_name:
        s.save()
    # Always return an HttpResponseRedirect after successfully dealing with POST data.
    # This prevents data from being posted twice if the User hits the back button.
    return HttpResponseRedirect('/strains/')
    ## save_new_strain ##


def add_contigs(request, strain_id):
    theStrain = Strain.objects.filter(pk = strain_id)[0]
    return render_to_response('strains/add_contigs.html',
                              { 'the_strain' : theStrain },
                              context_instance = RequestContext(request))
    ## add_contigs ##


def save_contigs(request, strain_id):
    theStrain = Strain.objects.filter(pk = strain_id)[0]
    contigFile = request.FILES['contig_file']
    contigs = []
    for seqRecord in SeqIO.parse(contigFile, "fasta"):
        c = Contig()
        c.strain = Strain.objects.filter(pk = strain_id)[0]
        c.name = seqRecord.id

        # The next line is a total kludge to the apparent limit of 1,000,000
        # (or so) characters in a models.TextField on MySQL. The nested calls
        # make a string of dna letters shrink by more than 50%. Since I don't
        # think there are any 2 megabase chromosomes in Sac cer this should
        # work.
        c.seq = base64.b64encode(bz2.compress(seqRecord.seq.tostring()))

        c.createdDate = datetime.datetime.now()
        c.modifiedDate = datetime.datetime.now()

        c.save()
    return HttpResponseRedirect('/strains/')
    ## save_new_contigs ##


def contigs(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    contigList = Contig.objects.filter(strain = theStrain).order_by('name')
    contigs = []
    for aContig in contigList:
        nFeatures = len(Feature.objects.filter(contig = aContig))
        nCDS = len(CDS.objects.filter(feature__contig = aContig))
        contigLen = len(bz2.decompress(base64.b64decode(aContig.seq)))
        contigs.append([aContig, contigLen, nFeatures, nCDS])

    return render_to_response('strains/contigs.html',
                              { 'the_strain'  : theStrain,
                                'contig_list' : contigs })
    ## contigs ##


def contig_detail(request, contig_id):
    theContig = get_object_or_404(Contig, pk = contig_id)
    theStrain = theContig.strain
    contigName = theContig.name

    # Don't forget to turn the string we saved in the db back into a DNA
    # sequence
    contigSeq = bz2.decompress(base64.b64decode(theContig.seq))

    contigSeqParts = []
    n = 0
    size = 60
    while 1:
        part = contigSeq[n * size : (n + 1) * size]
        if not part:
            break
        contigSeqParts.append(part)
        n += 1
    contigSeq = '\n'.join(contigSeqParts)

    featureCount = len(Feature.objects.filter(contig = theContig))

    return render_to_response('strains/contig_detail.html',
                              { 'the_contig'    : theContig,
                                'the_strain'    : theStrain,
                                'contig_seq'    : contigSeq,
                                'feature_count' : featureCount })
    ## contigDetail ##


def chop60(seq):
    seqParts = []
    n = 0
    size = 60
    while 1:
        part = seq[n * size : (n + 1) * size]
        if not part:
            break
        seqParts.append(part)
        n += 1
    return '\n'.join(seqParts)
    ## chop60 ##


def add_features(request, strain_id):
    theStrain = Strain.objects.filter(pk = strain_id)[0]
    return render_to_response('strains/add_features.html',
                              { 'the_strain' : theStrain },
                              context_instance = RequestContext(request))
    ## add_features ##
    

def save_features(request, strain_id):
    """
    seqid = line[0]
    source = line[1]
    type = line[2]
    start = line[3]
    end = line[4]
    score = line[5]
    strand = line[6] # "+", "-", or "."
    phase = line[7]
    attributes = line[8]
    """
    theStrain = Strain.objects.get(pk = strain_id)
    gffFile = request.FILES['gff_file']
    gff = gffFile.read()
    gff = gff.split('###')[0] # Throw away the sequence
    gff = [x.split('\t') for x in gff.splitlines() if x[0] != '#'] # Throw away the header comments. Now we're left with just the meat of the file

    contigMap = {}
    for seqid, source, featureType, start, end, score, strand, phase, attributes in gff:
        attributeParts = attributes.split(';')
        attributeParts = [x.split('=') for x in attributeParts]
        attributeParts = [(x[0], x[1].split(',')) for x in attributeParts]
        attributeParts = [(x[0], [urllib.unquote(y) for y in x[1]]) for x in attributeParts]

        attributeDict = {}
        for key, value in attributeParts:
            attributeDict[key] = value

        if featureType == 'contig':
            # We need to add this to the contigMap
            contigName = attributeDict['dbxref'][0].split(':')[-1]
            contigMap[seqid] = contigName
        else: # This is an actual feature line. It is assumed that we have already gone through all the contig lines
            theContig = get_object_or_404(Contig, name=contigMap[seqid] ) # Get the Contig we're going to point to
            feature = Feature()
            feature.contig = theContig
            try:
                feature.feature_id = attributeDict['ID'][0]
            except KeyError:
                pass
            else:
                # This one has a name that might be found in the Reference table
                if feature.feature_id.find(theStrain.name) != -1:
                    # Yup, it's one we need to link to the Reference table
                    referenceName = feature.feature_id.split("_")[0]
                    feature.reference = Reference.objects.get(feature_name = referenceName)
            if 'Parent' in attributeDict:
                parent = get_object_or_404(Feature, feature_id = attributeDict['Parent'][0], contig = theContig)
                feature.parent = parent
            feature.feature_type = featureType
            feature.start_coord = int(start)
            feature.stop_coord = int(end)
            feature.strand = strand

            feature.createdDate = datetime.datetime.now()
            feature.modifiedDate = datetime.datetime.now()

            feature.save()
    return HttpResponseRedirect('/strains/')
    ## save_features ##


def features_by_contig(request, contig_id):
    theContig = get_object_or_404(Contig, pk = contig_id)
    theStrain = theContig.strain
    features = Feature.objects.filter(contig = theContig).order_by('start_coord')
    features = [x for x in features if x.feature_id and not x.parent]
    for i, feature in enumerate(features):
        features[i] = [feature, get_child_features(feature)]

    return render_to_response('strains/features_by_contig.html',
                              { 'the_contig' : theContig,
                                'the_strain' : theStrain,
                                'feature_list' : features })
    ## features_by_contig ##


def get_child_features(feature):
    """
    get_child_features(feature):

    Given a Feature instance, finds it's children and returns them as a list of
    childList sublists. It calls itself recursively to fill in the child's
    children. When a feature has no children it returns an empty list.
    """
    children = []
    childFeatures = Feature.objects.filter(parent = feature).order_by('start_coord')
    if childFeatures:
        for childFeature in childFeatures:
            children.append([childFeature, get_child_features(childFeature)])
    return children
    ## get_child_features ##


def features_by_strain(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    contigList = Contig.objects.filter(strain = theStrain).order_by('name')
    featureList = []
    for aContig in contigList:
        contigFeatures = Feature.objects.filter(contig = aContig).order_by('start_coord')
        if len(contigFeatures) > 0:
            featureList.append([aContig, contigFeatures])
    return render_to_response('strains/features_by_strain.html',
                              { 'the_strain'   : theStrain,
                                'feature_list' : featureList })
    ## features_by_strain ##


def feature_detail(request, feature_id):
    theFeature = get_object_or_404(Feature, pk = feature_id)
    theContig = theFeature.contig

    contigSeq = bz2.decompress(base64.b64decode(theContig.seq))
    featureSeq = contigSeq[theFeature.start_coord - 1: theFeature.stop_coord]

    if theFeature.strand == '-':
        bpSeq = Seq(featureSeq, IUPAC.unambiguous_dna)
        featureSeq = bpSeq.reverse_complement().tostring()
    featureSeq = chop60(featureSeq)

    if not theFeature.parent: # This is a top_level feature
        featureName = theFeature.feature_id
    else:
        featureType = theFeature.feature_type
        aFeature = theFeature
        while aFeature.parent:
            aFeature = aFeature.parent
        featureName = '%s_%s' % (aFeature.feature_id, featureType)

    return render_to_response('strains/feature_detail.html',
                              { 'feature_seq'  : featureSeq,
                                'feature_name' : featureName })
    ## feature_detail ##


def add_CDS(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    return render_to_response('strains/add_CDS.html',
                              { 'the_strain' : theStrain },
                              context_instance = RequestContext(request))                              
    ## add_CDS ##


def add_protein(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    return render_to_response('strains/add_protein.html',
                              { 'the_strain' : theStrain },
                              context_instance = RequestContext(request))                              
    ## add_protein ##


def save_CDS(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    CDS_file = request.FILES['cds_file']
    for seqRecord in SeqIO.parse(CDS_file, "fasta"):
        cds = CDS()
        cds.name = seqRecord.id
        cds.seq = seqRecord.seq.tostring()
        cds.feature = Feature.objects.get(feature_id = cds.name)
        cds.createdDate = datetime.datetime.now()
        cds.modifiedDate = datetime.datetime.now()
        cds.save()
    return HttpResponseRedirect('/strains/')
    ## save_CDS ##


def save_protein(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    protein_file = request.FILES['protein_file']
    for seqRecord in SeqIO.parse(protein_file, "fasta"):
        protein = Protein()
        protein.name = seqRecord.id
        protein.seq = seqRecord.seq.tostring()
        protein.cds = CDS.objects.get(name = protein.name)
        protein.createdDate = datetime.datetime.now()
        protein.modifiedDate = datetime.datetime.now()
        protein.save()
    return HttpResponseRedirect('/strains/')
    ## save_protein ##


def references(request):
    referenceList = Reference.objects.filter(feature_type = 'ORF').order_by('standard_name')
    references = []
    for ref in referenceList:
        references.append([ref, len(CDS.objects.filter(feature__reference = ref))])
    return render_to_response('strains/references.html',
                              { 'reference_list' : references })
    ## references ##


def gene_detail(request, reference_id):
    theReference = get_object_or_404(Reference, pk = reference_id)
    theFeatures = Feature.objects.filter(reference = theReference)
    cdsGroup = CDS.objects.filter(feature__reference = theReference)
    proteinGroup = Protein.objects.filter(cds__feature__reference = theReference)

    cdsList = []
    for cds in cdsGroup:
        cdsList.append([cds, chop60(cds.seq)])
    proteinList = []
    for pep in proteinGroup:
        proteinList.append([pep, chop60(pep.seq)])
    return render_to_response('strains/gene_detail.html',
                              { 'the_reference' : theReference,
                                'feature_list'  : theFeatures,
                                'cds_list'      : cdsList,
                                'protein_list'  : proteinList })
    ## gene_detail ##


def cds_clustal(request, reference_id):
    return run_clustal(request, reference_id, 'cds')
    ## cds_clustal ##


def protein_clustal(request, reference_id):
    return run_clustal(request, reference_id, 'protein')
    ## cds_clustal ##


def run_clustal(request, reference_id, cdsOrPep):
    theReference = get_object_or_404(Reference, pk = reference_id)
    if cdsOrPep == 'cds':
        objList = CDS.objects.filter(feature__reference = theReference)
    else:
        objList = Protein.objects.filter(cds__feature__reference = theReference)

    handle, path = tempfile.mkstemp(suffix = '.fasta', prefix = 'strains_', text = True)
    outf = open(path, 'w')
    for obj in objList:
        outf.write(">%s\n%s\n" % (obj.name, chop60(obj.seq)))
    outf.close()

    command = ClustalwCommandline("clustalw2", infile = path)
    stdout, stderr = command()
    os.unlink(path)
    inf = open(path.replace(".fasta", ".aln"))
    result = inf.read()
    inf.close()
    if not result:
        result = 'No clustal alignment'
    os.unlink(path.replace(".fasta", ".aln"))
    try:
        os.unlink(path.replace(".fasta", ".dnd"))
    except OSError: # There was no .dnd file created. That's OK, go on
        pass
    return render_to_response('strains/alignment.html',
                              { 'the_reference' : theReference,
                                'result'        : result,
                                'cdsOrPep'      : cdsOrPep })
    ## run_clustal ##
    

def strain_cds(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    cdsList = CDS.objects.filter(feature__contig__strain = theStrain)
    thingList = []
    for cds in cdsList:
        thingList.append([cds.name, chop60(cds.seq)])
    return render_to_response('strains/fasta.html',
                              { 'things' : thingList })
    ## strain_cds ##


def strain_protein(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    pepList = Protein.objects.filter(cds__feature__contig__strain = theStrain)
    thingList = []
    for pep in pepList:
        thingList.append([pep.name, chop60(pep.seq)])
    return render_to_response('strains/fasta.html',
                              { 'things' : thingList })
    ## strain_protein ##


def strain_contigs(request, strain_id):
    theStrain = get_object_or_404(Strain, pk = strain_id)
    contigList = Contig.objects.filter(strain = theStrain)
    thingList = []
    for contig in contigList:
        thingList.append([contig.name, chop60(bz2.decompress(base64.b64decode(contig.seq)))])
    return render_to_response('strains/fasta.html',
                              { 'things' : thingList })
    ## strain_protein ##
