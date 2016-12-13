#!/nfs/goldstein/software/python2.7.7/bin/python

# append_vep.py

import re
from optparse import OptionParser
import subprocess
from tempfile import mkdtemp, mkstemp
import csv
import os
from collections import defaultdict
from fnmatch import fnmatch
import shutil
import pdb

parser = OptionParser()
parser.add_option("-i", "--atav_file",
                  dest="atav",
                  help="The CSV file generated by ATAV trio analysis.\
                       \nWiki: https://redmine.igm.cumc.columbia.edu/projects/bioinfo_tools/wiki/Diagnostic_Analysis_Framework")

parser.add_option("--append_file",
                  dest="output_",
                  help="The file path to to write the appended CSV to.")

parser.add_option("--output_file",
                  dest="vep",
                  help="OPTION: The position to write the VEP output file. If not given, the raw VEP output is not written.")

parser.add_option("--force_overwrite",
		action="store_true", dest="forceoverwrite", default=False,
		help = "overwrite the output file specified by --output_file")

parser.add_option("--force_append",
		action="store_true", dest="forceappend", default=False,
		help = "overwrite the appended file specified by --appended_file")

class Variant(object):
    eff_key = ['intergenic_variant', 'feature_truncation', 'regulatory_region_variant', 'feature_elongation',
               'regulatory_region_amplification', 'regulatory_region_ablation', 'TF_binding_site_variant',
               'TFBS_amplification', 'TFBS_ablation', 'downstream_gene_variant', 'upstream_gene_variant',
               'nc_transcript_variant', 'NMD_transcript_variant', 'intron_variant', 'non_coding_exon_variant',
               '3_prime_UTR_variant', '5_prime_UTR_variant', 'mature_miRNA_variant', 'coding_sequence_variant',
               'synonymous_variant', 'stop_retained_variant', 'incomplete_terminal_codon_variant',
               'splice_region_variant', 'protein_altering_variant', 'missense_variant', 'inframe_deletion',
               'inframe_insertion', 'transcript_amplification', 'initiator_codon_variant',
               'start_lost', 'stop_lost', 'frameshift_variant',
               'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant', 'transcript_ablation']
    ties = []
    def __init__(self, ln):
        ln = ln.strip().split('\t')
        self.line=ln

        self.id_ =  ln[0]

        self.gene = ln[3]
        self.feature = ln[4]
        self.feature_type = ln[5]
        self.consequences = ln[6].split(',')

        # process whitespace-free string at the end of variant data
        jargon = ln[-1].strip().split(';')
        self.polyphen = self.find_polyphen(jargon)
        self.hgvsc, self.hgvsp, self.exon = self.find_hgvs(jargon)
        self.ccds = self.find_ccds(jargon)
        self.canonical = False
        if 'canonical=' in str(jargon).lower():
            self.canonical = True

        h = self.consequences[0]
        for con in self.consequences[1:]:
            try:
                if Variant.eff_key.index(con) > Variant.eff_key.index(h):
                    h = con
            except ValueError:
                if not options.mute:
                    raise Warning ('Encountered un-indexed consequence: %s' % (self.con) )
        self.sev_con = Variant.eff_key.index(h)
        self.con = h

    def find_polyphen(self, lst):
        """returns PolyPhen score of Variant
        if none found, returns 0"""

        for field in lst:
            if 'polyphen' in field.lower():
                score = re.sub(r'[^\d.]+', '', field)
                return float(score)
        return 0

    def find_hgvs(self, lst):
        """returns HGVSc, HGVSp in tuple form
        if one score is found, returns a None in place

            ex. if HGVSp not found
            returns (HGVSc, None)"""
        hgvsc, hgvsp, exon = None, None, None
        for field in lst:
            if 'hgvsc=' == field.lower()[:6]:
                hgvsc= field[6:]
            if 'hgvsp=' == field.lower()[:6]:
                hgvsp= field[6:]
            if 'exon=' == field.lower()[:5]:
                exon= field[5:]
        return hgvsc, hgvsp, exon

    def find_ccds(self, lst):
        """returns a boolean representing whether the Variant has a CCDS
        """
        for field in lst:
            if 'ccds=' == field.lower()[:5]:
                return True
        return False

    def setfreq(self, n, d):
        self.freq = "%s/%s" % ( n , d )

    def compare (var1, var2):
        """
        Compares Variants and returns the most deleterious.

        Parameters
        ----------
        var1, var2: the Variants to compare

        Returns
        -------
        Variant: returns an instance of the Variant class determined
        the following prioritization of criteria
            1. most damaging consequence
            2. highest PolyPhen score
            3. the canonical variant, if one variant is canonical and the var2 is not
            4. the transcript that has a CCDS ID, if the other does not
            5. var1
        """

        if var1.sev_con > var2.sev_con:
            return var1
        elif var1.sev_con == var2.sev_con:
            if var1.polyphen > var2.polyphen:
                return var1
            elif var1.polyphen == var2.polyphen:
                if var1.canonical and not var2.canonical:
                    return var1
                elif var1.canonical is var2.canonical:
                    if var1.ccds and not var2.ccds:
                        return var1
                    elif var1.ccds and var2.ccds:
                        if var1.feature != var2.feature and var1.canonical and var2.canonical:
			    Variant.ties.append(var1.feature)
			    Variant.ties.append(var2.feature)
                            return var1
        return var2

    def same_pos(self, other):
	if other.chrom != self.chrom:
	    return False

	elif len(self.position) == 1 and len(other.position) == 1:
	    if self.position == other.position:
		return True
	    return False
	else:
	    return not set(self.position).isdisjoint(set(other.position))


class csvRow(list):
    """A list which can be written directly to a csv file"""
    def __init__(self, *args):
        list.__init__(self, *args)
        if self[-1] is '':
            del(self[-1])
    def append(self, p_object):
        if p_object is None:
            p_object= 'NA'
        super(csvRow, self).append(p_object)
    def __str__(self):
        s=''
        for c in self:
            s+=c+','
        return s+'\n'


def doVEP(fl):
    # first get vep input from atav output file
    i = open(fl, 'rb')

    reader = csv.reader(i, delimiter=',')

    # find what column the variant IDs are in
    col = []
    header = reader.next()
    for cell in header:
        if fnmatch(cell.lower(),"*variant*id*"):
            col.append(header.index(cell))
    if col == None:
        raise NameError('No Variant ID Column Found in input.')

    variantids = []

    for row in reader:
	variantids += [row[c] for c in col]
    i.close()

    variantids = list(set(variantids))

    if options.vep is None:
	out, out_path = mkstemp()
    else:
        out_path = options.vep

    inp, inp_path = mkstemp()

    lin = [v.split('-') for v in variantids]
    try:
        for v in variantids:
            var = v.split('-')
            out_form = var[0] + '\t' + var[1] + '\t' + v + '\t' + var[2] + '\t' + var[3] + '\n'
            os.write(inp, out_form)

	vep_call = ["perl", "/nfs/goldstein/goldsteinlab/software/variant_effect_predictor_74/variant_effect_predictor.pl",
                                 "-i", inp_path, "--format", "vcf", "--dir", "/nfs/goldstein/goldsteinlab/software/variant_effect_predictor_74",
                                 "-o", out_path, "--ccds", " --protein", "--domains", "--polyphen=s", "--sift=b", "--symbol",
                                 "--numbers",  "--hgvs", "--canonical",  "--cache", "--offline", "--fasta",
                                 "/nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa"]
	if options.forceoverwrite or not options.vep:
	    vep_call.append('--force_overwrite')
	if len(variantids) > 0 :
		try:
		    p = subprocess.check_output(vep_call)
		except subprocess.CalledProcessError as e:
		    raise EnvironmentError("could not run perl VEP command\n{}".format(" ".join(vep_call)))

        f = open(out_path)
        lines = f.readlines()
        f.close()
    finally:
	os.remove(inp_path)
	if options.vep is None:
	    os.remove(out_path)

    variants = defaultdict(list)
    for l in lines:
	if l[0] == '#':
	    continue
	else:
	    lineVariant = Variant(l)
        variants[lineVariant.id_].append(lineVariant)

    for i in range(len(variantids)):
	a = variantids[i].split('-')[:-2]
	variantids[i] = ( a[0] , int(a[1]) )

    return variants, variantids


def readVEP(variants, atavids):
    """Phase 1: parses through the VEP file and returns the three output columns
    for each variant in lists (see returns for format)

    Parameters
    ----------
    arg: the path to the VEP output file as a string

    Returns
    -------
    (variants, col3): where hgvs and col3 are both lists
        variants is a list Variant instances representing the most damaging variants
        for each Variant ID in the VEP file

        col3 is n/d for each variant (hence it will always be half the size of hgvs)
            where n is the number of transcripts of the top gene that have the most damaging consequence
            and d is the number the number of transcripts in the gene"""

    highest = {}

    for varid, variant_cons in variants.iteritems():
        h = variant_cons[0]
        for variant in variant_cons:
            h = Variant.compare(h, variant)

        highest[ h.id_ ] = h
        n = sum(1 for v in variant_cons if Variant.eff_key[h.sev_con] in v.consequences
                and h.gene == v.gene and v.feature_type == 'Transcript' and v.ccds)
        d = sum(1 for v in variant_cons if h.gene == v.gene and v.feature_type == 'Transcript' and v.ccds)
        h.setfreq(n , d)

    return highest


def makeOutput(arg, highest):
    reader = open(arg, 'r')

    # first make header
    header = reader.next().strip().split(',')
    id_col = [col_num for col_num,cell in enumerate(header) if fnmatch(cell.lower(),"*variant*id*")]
    header = csvRow(header)
    if header[-1] == '': del header[-1]
    num = ""
    if len(id_col) > 1: num = " (#1)"

    cols = 4 # number of columns of output
    header.append('VEP Function{}'.format(num))
    header.append('HGVSc{}'.format(num))
    header.append('HGVSp{}'.format(num))
    header.append('Exon{}'.format(num))
    header.append('#Transcripts{}'.format(num))

    if len(id_col) == 2:
        num = " (#2)"
        header.append('VEP Function{}'.format(num))
        header.append('HGVSc{}'.format(num))
        header.append('HGVSp{}'.format(num))
        header.append('Exon{}'.format(num))
        header.append('#Transcripts{}'.format(num))


    rows = [csvRow(header)]
    for row in reader:
        row = row.strip().split(',')
        row = csvRow(row)
        transcript = highest[ row[id_col[0]].strip('"') ]
        row.append(transcript.con)
        row.append(transcript.hgvsc)
        row.append(transcript.hgvsp)
        row.append(transcript.exon)
        row.append(transcript.freq)
        if len(id_col) == 2:
            transcript = highest[ row[id_col[1]].strip('"') ]
            row.append(transcript.con)
            row.append(transcript.hgvsc)
            row.append(transcript.hgvsp)
            row.append(transcript.exon)
            row.append(transcript.freq)
        rows.append(row)
    reader.close()

    return rows

def writeOutput(arg, rows):
    o = open(arg, 'w+')
    for row in rows:
        o.write(str(row))
    o.close()

def append_vep():
    if not options.atav:
        parser.error('ATAV CSV input file not given. Specify with -i')
    if not options.output_:
        parser.error('Output file not given. Specify with --appended_file')

    if os.path.exists(options.output_) and not options.forceappend:
	raise IOError( "Error: Append file " + options.output_ + " already exists. Specify a different output file with --append_file or overwrite existing file with --force_append.")

    print " - Running VEP..."
    vep, atavids = doVEP(options.atav)
    print " - Reading VEP..."
    del_var = readVEP(vep, atavids)
    print " - Making output..."
    rows = makeOutput(options.atav, del_var)
    print " - Writing output..."
    writeOutput(options.output_, rows)
    if options.vep:
	print " - Wrote stats to " + options.vep + "_summary.html"
    print " - Done!"


options, args = parser.parse_args()
append_vep()

