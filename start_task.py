import sys
sys.path.append('/opt/boglab')
sys.path.append('/opt/boglab/genome_scoring')
from datetime import datetime
from optparse import OptionParser
from genome_scoring import TASK_MODULE
import re
import urllib

# import arguments and options
usage = 'usage: %prog [options]'
parser = OptionParser(usage=usage)
parser.add_option('-f', '--fasta', dest='fasta', type='string', default='NA', help='Legacy.')
parser.add_option('-r', '--rvds', dest='rvds', type = 'string', default='NA', help='RVD sequence seperated by spaces.')
parser.add_option('-o', '--organism', dest='organism', type = 'string', default='NA', help='Name of organism for the genome to be searched.')
parser.add_option('-p', '--outpath', dest='outpath', type='string', default = 'NA', help='Optional full path for output file; if set --job, --outdir and --outfile are ignored.')
parser.add_option('-y', '--offtarget', dest='offtarget', action = 'store_true', default = False, help='Legacy')
parser.add_option('-z', '--nodeid', dest='nodeid', type='int', default = '-1', help='Optional node id if this script was called from Drupal.')
(options, args) = parser.parse_args()

def logger(message):
    print "[%s] %s" % (datetime.now().ctime(), message)


RVD_re = re.compile(r'^(?:(NI|NN|NG|HD|NS)\s*){12,20}$', re.IGNORECASE | re.MULTILINE)
if not RVD_re.match(options.rvds):
    urllib.urlopen("https://boglab.plp.iastate.edu/talent/jobcomplete/" + str(options.nodeid) + "/1")
    logger("RVD sequence is not in the correct format.  Enter between 12 and 20 RVDs using the standard single letter amino acid abbreviations. Currently supported RVDs for genomes are HD, NN, NI, NG, NS.")
    sys.exit(1)

valid_organisms = ['drosophila_melanogaster', 'arabidopsis_thaliana', 'mus_musculus', 'oryza_sativa', 'caenorhabditis_elegans']

if options.organism not in valid_organisms:
    urllib.urlopen("https://boglab.plp.iastate.edu/talent/jobcomplete/" + str(options.nodeid) + "/1")
    logger("Invalid organism specified.")
    sys.exit(1)
    
_temp = __import__(TASK_MODULE % options.organism, globals(), locals(), ['TalentTask'], -1)
TalentTask = _temp.TalentTask

logger("Working...")

TalentTask.apply_async(args=(options.rvds, options.outpath, options.nodeid), queue=options.organism)
