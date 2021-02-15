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

import csv
import traceback
import pytransit.transit_tools as transit_tools
from pytransit.convert import base

class InvalidArgumentException(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        super(InvalidArgumentException, self).__init__(message)

############# Description ##################

short_name = "gff_to_prot"
long_name = "GFF to Prot_table converter"
description = "Convert a GFF file to Prot_table format"
label = "GFF3 to Prot_table"

############# Analysis Method ##############

class GffProtConverter(base.TransitConvert):
    def __init__(self):
        base.TransitConvert.__init__(self, short_name, long_name, description, label, GffProtMethod, GffProtGUI)


################# GUI ##################
class GffProtGUI(base.ConvertGUI):

    def __init__(self):
        base.ConvertGUI.__init__(self)

########## METHOD #######################

class GffProtMethod(base.ConvertMethod):
    """
    GffProtMethod
    """
    def __init__(self,
                annotation_path,
                output, wxobj=None):
        self.short_name = short_name
        self.long_name = long_name
        self.description = description
        self.label = label
        self.output = output
        self.annotation_path = annotation_path
        base.ConvertMethod.__init__(self, short_name, long_name, description, label, annotation_path, output, wxobj=wxobj)

    @classmethod
    def fromGUI(self, wxobj):
        """ """
        #Get Annotation file
        annotationPath = wxobj.annotation
        if not transit_tools.validate_annotation(annotationPath):
            return None

        #Get output path
        defaultFileName = "{0}.prot_table".format(os.path.splitext(os.path.basename(annotationPath))[0])
        defaultDir = os.getcwd()
        output_path = wxobj.SaveFile(defaultDir, defaultFileName)
        if not output_path: return None
        output_file = open(output_path, "w")

        return self(annotationPath, output_file, wxobj)

    @classmethod
    def fromargs(self, rawargs):
        (args, kwargs) = transit_tools.cleanargs(rawargs)
        if (len(args) < 2):
            print("Error: Please specify Input and Output paths")
            print(self.usage_string())
            sys.exit(1)

        annotationPath = args[0]
        outpath = args[1]
        output_file = open(outpath, "w")

        return self(annotationPath, output_file)

    @classmethod
    def fromconsole(self):
        try:
            return self.fromargs(sys.argv[3:])
        except InvalidArgumentException as e:
            print("Error: %s" % str(e))
            print(self.usage_string())
        except IndexError as e:
            print("Error: %s" % str(e))
            print(self.usage_string())
        except TypeError as e:
            print("Error: %s" % str(e))
            traceback.print_exc()
            print(self.usage_string())
        except ValueError as e:
            print("Error: %s" % str(e))
            traceback.print_exc()
            print(self.usage_string())
        except Exception as e:
            print("Error: %s" % str(e))
            traceback.print_exc()
            print(self.usage_string())
        sys.exit()

    def get_description(self, line, parent):
        cols = line.split('\t')
        labels = {}
        print(line)
        print(len(cols))
        for pair in cols[8].split(";"):
            k, v = pair.split('=')
            labels[k] = v

        if (cols[2]) == "CDS" and labels["Parent"] == parent:
            #return labels.get("Note", '-')
            return labels.get("product", '-')
        return '-'

    def Run(self):
        gff_file = open(self.annotation_path)
        output_file = self.output
        writer = csv.writer(output_file, delimiter='\t')
        lines = gff_file.readlines()
        gff_file.close()
        self.transit_message("Converting annotation file from GFF3 format to prot_table format")

        for i, line in enumerate(lines):
            line = line.strip()
            if len(line)==0 or line.startswith('#'): continue
            cols = line.split('\t')
            if (len(cols) < 9): 
                sys.stderr.write(("Ignoring invalid row with entries: {0}\n".format(cols)))
                continue
            if (cols[2]) == "CDS": # if you also want tRNAs and rRNAs, modify here
                if "locus_tag" not in line: print("warning: skipping lines that do not contain 'locus_tag'"); continue
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6].strip()
                size = int(abs(end-start+1)/3) # includes stop codon
                labels = {}
                for pair in cols[8].split(";"):
                    k, v = pair.split('=')
                    labels[k.strip()] = v.strip()
                Rv = labels["locus_tag"].strip() # error out if not found
                gene = labels.get('gene', '') # or Name?
                if gene=="": gene = '-'
                desc = labels.get('product', '') 
                vals = [desc, start, end, strand, size, '-', '-', gene, Rv, '-']
                writer.writerow(vals)
        output_file.close()
        self.transit_message("Finished conversion")

    @classmethod
    def usage_string(self):
        return """python %s convert gff_to_prot_table <annotation in gff format> <output file>""" % (sys.argv[0])

if __name__ == "__main__":

    pass

