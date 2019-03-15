# Copyright 2015.
#   Michael A. DeJesus, Chaitra Ambadipudi, and  Thomas R. Ioerger.
#
#
#    This file is part of TRANSIT.
#
#    TRANSIT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#
#    TRANSIT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.


import sys
import os


try:
    import wx
    WX_VERSION = int(wx.version()[0])
    hasWx = True

except Exception as e:
    hasWx = False
    WX_VERSION = 0

if hasWx:
    import wx.xrc
    from wx.lib.buttons import GenBitmapTextButton
    from pubsub import pub
    import wx.adv

import math
import ntpath
import numpy
import scipy.optimize
import scipy.stats

import warnings

import pytransit
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools


if hasWx:
    class AssumeZerosDialog(wx.Dialog):

        def __init__(self, *args, **kw):

            self.ID_HIMAR1 = wx.NewId()
            self.ID_TN5 = wx.NewId()

            wx.Dialog.__init__(self, None, title="Dialog")

            self.ID_HIMAR1 = wx.NewId()
            self.ID_TN5 = wx.NewId()
    
            self.SetSize((500, 300))
            self.SetTitle("Warning:  Wig Files Do Not Include Empty Sites")

            mainSizer = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(mainSizer)
    
            warningText = """

One or more of your .wig files does not include any empty sites (i.e. sites with zero read-counts). The analysis methods in TRANSIT require knowing ALL possible insertion sites, even those without reads.
    
    Please indicate how you want to proceed:

    As Himar1: You will need to provide the DNA sequence (.fasta format) and TRANSIT will automatically determine empty TA sites.

    As Tn5: TRANSIT will assume all nucleotides are possible insertion sites. Those not included in the .wig file are assumed to be zero.
    """
            warningStaticBox = wx.StaticText(self, wx.ID_ANY, warningText, (-1,-1), (-1, -1), wx.ALL)
            warningStaticBox.Wrap(480)
            mainSizer.Add(warningStaticBox, flag=wx.CENTER, border=5)
    
            button_sizer = wx.BoxSizer(wx.HORIZONTAL)
            himar1Button = wx.Button(self, self.ID_HIMAR1, label='Proceed as Himar1')
            tn5Button = wx.Button(self, self.ID_TN5, label='Proceed as Tn5')
            cancelButton = wx.Button(self, wx.ID_CANCEL, label='Cancel')

    
            button_sizer.Add(himar1Button, flag=wx.LEFT, border=5)
            button_sizer.Add(tn5Button, flag=wx.LEFT, border=5)
            button_sizer.Add(cancelButton, flag=wx.LEFT, border=5)

            mainSizer.Add(button_sizer,
                flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)


            himar1Button.Bind(wx.EVT_BUTTON, self.OnClose)
            tn5Button.Bind(wx.EVT_BUTTON, self.OnClose)
            cancelButton.Bind(wx.EVT_BUTTON, self.OnClose)


        def OnClose(self, event):
    
            if self.IsModal():
                self.EndModal(event.EventObject.Id)
            else:
                self.Close()



def aton(aa):
    #TODO: Write docstring
    return(((aa-1)*3)+1)

def parseCoords(strand, aa_start, aa_end, start, end):
    #TODO: Write docstring
    if strand == "+":
        return((aton(aa_start) + start,  aton(aa_end) + start))
    # Coordinates are Reversed... to match with Trash FILE TA coordinates
    if strand == "-":
        return((end - aton(aa_end), end - aton(aa_start)))



def fetch_name(filepath):
    #TODO: Write docstring
    return os.path.splitext(ntpath.basename(filepath))[0]


def basename(filepath):
    #TODO: Write docstring
    return ntpath.basename(filepath)


def dirname(filepath):
    return os.path.dirname(os.path.abspath(filepath))


def cleanargs(rawargs):
    """Returns a list and a dictionary with positional and keyword arguments.

    -This function assumes flags must start with a "-" and and cannot be a 
        number (but can include them).
    
    -Flags should either be followed by the value they want to be associated 
        with (i.e. -p 5) or will be assigned a value of True in the dictionary.

    -The dictionary will map flags to the name given minus ONE "-" sign in
        front. If you use TWO minus signs in the flag name (i.e. --verbose), 
        the dictionary key will be the name with ONE minus sign in front 
        (i.e. {"-verbose":True}).
    

    Arguments:
        rawargs (list): List of positional/keyword arguments. As obtained from
                         sys.argv.

    Returns:
        list: List of positional arguments (i.e. arguments without flags),
                in order provided.
        dict: Dictionary mapping flag (key is flag minus the first "-") and
                their values.

    """
    args = []
    kwargs = {}
    count = 0
    # Loop through list of arguments
    while count < len(rawargs):
        # If the current argument starts with "-", then it's probably a flag
        if rawargs[count].startswith("-"):
            # Check if next argument is a number
            try:
                temp = float(rawargs[count+1])
                nextIsNumber = True
            except:
                nextIsNumber = False

            stillNotFinished = count + 1 < len(rawargs)
            if stillNotFinished:
                nextIsNotArgument = not rawargs[count+1].startswith("-")
                nextLooksLikeList = len(rawargs[count+1].split(" ")) > 1
            else:
                nextIsNotArgument = True
                nextLooksLikeList = False


            # If still things in list, and they look like arguments to a flag, add them to dict 
            if stillNotFinished and (nextIsNotArgument or nextLooksLikeList or nextIsNumber):
                kwargs[rawargs[count][1:]] = rawargs[count+1]
                count += 1
            # Else it's a flag but without arguments/values so assign it True
            else:
                kwargs[rawargs[count][1:]] = True
        # Else, it's probably a positional arguement without flags
        else:
            args.append(rawargs[count])
        count += 1
    return (args, kwargs)


def getTabTableData(path, colnames):
    #TODO: Write docstring
    row = 0
    data = []
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        tmp[-1] = tmp[-1].strip()
        rowdict = dict([(colnames[i], tmp[i]) for i in range(len(colnames))])
        data.append((row, rowdict))
        row+=1

    return data


def ShowMessage(MSG=""):
    #TODO: Write docstring
    wx.MessageBox(MSG, 'Info',
        wx.OK | wx.ICON_INFORMATION)

def ShowAskWarning(MSG=""):
    #TODO: Write docstring
    dial = wx.MessageDialog(None, MSG, 'Warning',
        wx.OK | wx.CANCEL | wx.ICON_EXCLAMATION)
    return dial.ShowModal()

def ShowError(MSG=""):
    #TODO: Write docstring
    dial = wx.MessageDialog(None, MSG, 'Error',
        wx.OK | wx.ICON_ERROR)
    dial.ShowModal()

def transit_message(msg="", prefix=""):
    #TODO: Write docstring
    if prefix:
        print prefix, msg
    else:
        print pytransit.prefix, msg

def transit_error(text):
    #TODO: Write docstring
    transit_message(text)
    try:
        ShowError(text)
    except:
        pass


def validate_annotation(annotation):
    #TODO: Write docstring
    if not annotation or not os.path.exists(annotation):
        transit_error("Error: No or Invalid annotation file selected!")
        return False
    return True

def validate_control_datasets(ctrldata):
    #TODO: Write docstring
    if len(ctrldata) == 0:
        transit_error("Error: No control datasets selected!")
        return False
    return True

def validate_both_datasets(ctrldata, expdata):
    #TODO: Write docstring
    if len(ctrldata) == 0 and len(expdata) == 0:
        transit_error("Error: No datasets selected!")
        return False
    elif len(ctrldata) == 0:
        transit_error("Error: No control datasets selected!")
        return False
    elif len(expdata) == 0:
        transit_error("Error: No experimental datasets selected!")
        return False
    else:
        return True


def validate_transposons_used(datasets, transposons, justWarn=True): 

    #TODO: Write docstring
    # Check if transposon type is okay.
    unknown = tnseq_tools.get_unknown_file_types(datasets, transposons)
    if unknown:
        if justWarn:
            answer = ShowAskWarning("Warning: Some of the selected datasets look like they were created using transposons that this method was not intended to work with: %s. Proceeding may lead to errors. Click OK to continue." % (",". join(unknown)))
            if answer == wx.ID_CANCEL:
                return False
            else:
                return True
        else:
            transit_error("Error: Some of the selected datasets look like they were created using transposons that this method was not intended to work with: %s." % (",". join(unknown)))
            return False

    return True



def validate_wig_format(wig_list, wxobj=None):
    # Check if the .wig files include zeros or not
    status = 0
    genome = ""
    includesZeros = tnseq_tools.check_wig_includes_zeros(wig_list)

    if sum(includesZeros) < len(includesZeros):
        # If console mode, just print a warning
        if not wxobj or not hasWx:
            warnings.warn("\nOne or more of your .wig files does not include any empty sites (i.e. sites with zero read-counts). Proceeding as if data was Tn5 (all other sites assumed to be zero)!\n")
            return (2, "")

        # Else check their decision
        dlg = AssumeZerosDialog()
        result = dlg.ShowModal()
        if result == dlg.ID_HIMAR1 and wxobj:
            status = 1
            # Get genome
            wc = u"Known Sequence Extensions (*.fna,*.fasta)|*.fna;*.fasta;|\nAll files (*.*)|*.*"
            gen_dlg = wx.FileDialog(wxobj, message="Save file as ...", defaultDir=os.getcwd(), defaultFile="", wildcard=wc, style=wx.FD_OPEN)
            if gen_dlg.ShowModal() == wx.ID_OK:
                genome = gen_dlg.GetPath()
            else:
                genome = ""

        elif result == dlg.ID_TN5:
            status = 2; genome = "" 
        else:
            status = 3; genome = ""
    return (status, genome)


def validate_filetypes(datasets, transposons, justWarn=True):
    validate_transposons_used(datasets, transposons, justWarn)


def get_pos_hash(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.
    
    Arguments:
        path (str): Path to annotation in .prot_table or GFF3 format.
    
    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    filename, file_extension = os.path.splitext(path)
    if file_extension.lower() in [".gff", ".gff3"]:
        return tnseq_tools.get_pos_hash_gff(path)
    else:
        return tnseq_tools.get_pos_hash_pt(path)
       

def get_extended_pos_hash(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.

    Arguments:
        path (str): Path to annotation in .prot_table or GFF3 format.

    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    filename, file_extension = os.path.splitext(path)
    if file_extension.lower() in [".gff", ".gff3"]:
        return tnseq_tools.get_extended_pos_hash_gff(path)
    else:
        return tnseq_tools.get_extended_pos_hash_pt(path)



def get_gene_info(path):
    """Returns a dictionary that maps gene id to gene information.
    
    Arguments:
        path (str): Path to annotation in .prot_table or GFF3 format.
    
    Returns:
        dict: Dictionary of gene id to tuple of information:
            - name
            - description
            - start coordinate
            - end coordinate
            - strand
            
    """
    filename, file_extension = os.path.splitext(path)
    if file_extension.lower() in [".gff", ".gff3"]:
        return tnseq_tools.get_gene_info_gff(path)
    else:
        return tnseq_tools.get_gene_info_pt(path)


def convertToIGV(self, dataset_list, annotationPath, path, normchoice=None):

    if not normchoice:
        normchoice = "nonorm"

    (fulldata, position) = tnseq_tools.get_data(dataset_list)
    (fulldata, factors) = norm_tools.normalize_data(fulldata, normchoice, dataset_list, annotationPath)
    position = position.astype(int)

    output = open(path, "w")
    output.write("#Converted to IGV with TRANSIT.\n")
    if normchoice != "nonorm":
        output.write("#Reads normalized using '%s'\n" % normchoice)

    output.write("#Files:\n#%s\n" % "\n#".join(dataset_list))
    output.write("#Chromosome\tStart\tEnd\tFeature\t%s\tTAs\n" % ("\t".join([transit_tools.fetch_name(D) for D in dataset_list])))
    chrom = transit_tools.fetch_name(annotationPath)

    for i,pos in enumerate(position):
        output.write("%s\t%s\t%s\tTA%s\t%s\t1\n" % (chrom, position[i], position[i]+1, position[i], "\t".join(["%1.1f" % fulldata[j][i] for j in range(len(fulldata))])))
    output.close()




def convertToCombinedWig(dataset_list, annotationPath, outputPath, normchoice="nonorm"):
    """Normalizes the input datasets and outputs the result in CombinedWig format.
    
    Arguments:
        dataset_list (list): List of paths to datasets in .wig format
        annotationPath (str): Path to annotation in .prot_table or GFF3 format.
        outputPath (str): Desired output path.
        normchoice (str): Choice for normalization method.
            
    """


    (fulldata, position) = tnseq_tools.get_data(dataset_list)
    (fulldata, factors) = norm_tools.normalize_data(fulldata, normchoice, dataset_list, annotationPath)
    position = position.astype(int)

    hash = get_pos_hash(annotationPath)
    rv2info = get_gene_info(annotationPath)

    output = open(outputPath, "w")
    output.write("#Converted to CombinedWig with TRANSIT.\n")
    if normchoice != "nonorm":
        output.write("#Reads normalized using '%s'\n" % normchoice)
        if type(factors[0]) == type(0.0):
            output.write("#Normalization Factors: %s\n" % "\t".join(["%s" % f for f in factors.flatten()]))
        else:
            output.write("#Normalization Factors: %s\n" % " ".join([",".join(["%s" % bx for bx in b]) for b in factors]))

    (K,N) = fulldata.shape
    output.write("#Files:\n")
    for f in dataset_list:
        output.write("#%s\n" % f)

    for i,pos in enumerate(position):
        #output.write("%-10d %s  %s\n" % (position[i], "".join(["%7.1f" % c for c in fulldata[:,i]]),",".join(["%s (%s)" % (orf,rv2info.get(orf,["-"])[0]) for orf in hash.get(position[i], [])])   ))
        output.write("%d\t%s\t%s\n" % (position[i], "\t".join(["%1.1f" % c for c in fulldata[:,i]]),",".join(["%s (%s)" % (orf,rv2info.get(orf,["-"])[0]) for orf in hash.get(position[i], [])])   ))
    output.close()



def convertToGeneCountSummary(dataset_list, annotationPath, outputPath, normchoice="nonorm"):
    """Normalizes the input datasets and outputs the result in CombinedWig format.
    
    Arguments:
        dataset_list (list): List of paths to datasets in .wig format
        annotationPath (str): Path to annotation in .prot_table or GFF3 format.
        outputPath (str): Desired output path.
        normchoice (str): Choice for normalization method.
            
    """

    (fulldata, position) = tnseq_tools.get_data(dataset_list) 
    (fulldata, factors) = norm_tools.normalize_data(fulldata, normchoice, dataset_list, annotationPath)
    output = open(outputPath, "w")
    output.write("#Summarized to Mean Gene Counts with TRANSIT.\n")
    if normchoice != "nonorm":
        output.write("#Reads normalized using '%s'\n" % normchoice)
        if type(factors[0]) == type(0.0):
            output.write("#Normalization Factors: %s\n" % "\t".join(["%s" % f for f in factors.flatten()]))
        else:
            output.write("#Normalization Factors: %s\n" % " ".join([",".join(["%s" % bx for bx in b]) for b in factors]))

    (K,N) = fulldata.shape
    output.write("#Files:\n")
    for f in dataset_list:
        output.write("#%s\n" % f)

    # Get Gene objects
    G = tnseq_tools.Genes(dataset_list, annotationPath, norm=normchoice)

    dataset_header = "\t".join([os.path.basename(D) for D in dataset_list])
    output.write("#Orf\tName\tNumber of TA sites\t%s\n" % dataset_header)
    for i,gene in enumerate(G):
        if gene.n > 0:
            data_str = "\t".join(["%1.2f" % (M) for M in numpy.mean(gene.reads, 1)])
        else:
            data_str = "\t".join(["%1.2f" % (Z) for Z in numpy.zeros(K)])
        output.write("%s\t%s\t%s\t%s\n" % (gene.orf, gene.name, gene.n, data_str))
    output.close()




def get_validated_data(wig_list, wxobj=None):
    """ Returns a tuple of (data, position) containing a matrix of raw read-counts
        , and list of coordinates. 

    Arguments:
        wig_list (list): List of paths to wig files.
        wxobj (object): wxPython GUI object for warnings

    Returns:
        tuple: Two lists containing data and positions of the wig files given.

    :Example:

        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> (data, position) = tnseq_tools.get_validated_data(["data/glycerol_H37Rv_rep1.wig", "data/glycerol_H37Rv_rep2.wig"])
        >>> print data
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])

    .. seealso:: :class:`get_file_types` :class:`combine_replicates` :class:`get_data_zero_fill` :class:`pytransit.norm_tools.normalize_data`
    """

    (status, genome) = validate_wig_format(wig_list, wxobj=wxobj)

    # Regular file with empty sites
    if status == 0:
        return tnseq_tools.get_data(wig_list)    
    # No empty sites, decided to proceed as Himar1
    elif status == 1:
        return tnseq_tools.get_data_w_genome(wig_list, genome)
    # No empty sites, decided to proceed as Tn5
    elif status == 2:
        return tnseq_tools.get_data_zero_fill(wig_list)
    # Didn't choose either.... what!?
    else:
        return tnseq_tools.get_data([])

