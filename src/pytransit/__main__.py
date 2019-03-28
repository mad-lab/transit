
import sys

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


import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.analysis
import pytransit.export
import pytransit.convert

method_wrap_width = 250
methods = pytransit.analysis.methods
export_methods = pytransit.export.methods
convert_methods = pytransit.convert.methods
all_methods = {}
all_methods.update(methods)

#all_methods.update(export_methods)

wildcard = "Python source (*.py)|*.py|" \
            "All files (*.*)|*.*"
transit_prefix = "[TRANSIT]"



def run_main():
    (args, kwargs) = transit_tools.cleanargs(sys.argv[1:])
    main(*args, **kwargs)

def main(*args, **kwargs):
    #If no arguments, show GUI:
    DEBUG = "--debug" in sys.argv
    if DEBUG:
        sys.argv.remove("--debug")
        kwargs.pop("-debug")

    if (not args and ('v' in kwargs or '-version' in kwargs)):
        print "Version: {0}".format(pytransit.__version__)
        sys.exit(0)
    if (not args and ('h' in kwargs or '-help' in kwargs)):
        print "For commandline mode, please use one of the known methods (or see documentation to add a new one):"
        print("Analysis methods: ")
        for m in all_methods:
            ## TODO :: Move normalize to separate subcommand?
            if (m == "normalize"): continue
            print "\t - %s" % m
        print("Other functions: ")
        print("\t - normalize")
        print("\t - convert")
        print("\t - export")
        print "Usage: python %s <method>" % sys.argv[0]
        sys.exit(0)

    # Check if running in GUI Mode
    if not (args or kwargs) and hasWx:

        import matplotlib
        matplotlib.use("WXAgg")
        import matplotlib.pyplot
        import pytransit.transit_gui as transit_gui
        transit_tools.transit_message("Running in GUI Mode")
        app = wx.App(False)

        #create an object of CalcFrame
        frame = transit_gui.TnSeekFrame(None, DEBUG)
        #show the frame
        frame.Show(True)
        frame.Maximize(True)

        #start the applications
        app.MainLoop()

    # Tried GUI mode but has no wxPython
    elif not (args or kwargs) and not hasWx:
        print "Please install wxPython to run in GUI Mode."
        print "To run in Console Mode please follow these instructions:"
        print ""
        print "Usage: python %s <method>" % sys.argv[0]
        print "List of known methods:"
        for m in methods:
            print "\t - %s" % m
    # Running in Console mode
    else:
        import matplotlib
        matplotlib.use("Agg")
        method_name = args[0]
        if method_name not in all_methods:
            if method_name.lower() == "export":
                export_method_name = ""
                if len(args) > 1:
                    export_method_name = args[1]

                if export_method_name not in export_methods:
                    print "Error: Need to specify the export method."
                    print "Please use one of the known methods (or see documentation to add a new one):"
                    for m in export_methods:
                        print "\t - %s" % m
                    print "Usage: python %s export <method>" % sys.argv[0]
                else:
                    methodobj = export_methods[export_method_name].method.fromconsole()
                    methodobj.Run()
            elif method_name.lower() == "convert":
                convert_method_name = ""
                if len(args) > 1:
                    convert_method_name = args[1]

                if convert_method_name not in convert_methods:
                    print "Error: Need to specify the convert method."
                    print "Please use one of the known methods (or see documentation to add a new one):"
                    for m in convert_methods:
                        print "\t - %s" % m
                    print "Usage: python %s convert <method>" % sys.argv[0]
                else:
                    methodobj = convert_methods[convert_method_name].method.fromconsole()
                    methodobj.Run()
            else:
                print "Error: The '%s' method is unknown." % method_name
                print "Please use one of the known methods (or see documentation to add a new one):"
                for m in all_methods:
                    print "\t - %s" % m
                print "Usage: python %s <method>" % sys.argv[0]
        else:

            methodobj = all_methods[method_name].method.fromconsole()
            methodobj.Run()



if __name__ == "__main__":
    main()
    sys.exit(1)

