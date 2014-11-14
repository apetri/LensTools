from ..simulations import Gadget2Snapshot
import sys
import argparse

def main(args=None):

	parser = argparse.ArgumentParser()
	parser.add_argument("filename",nargs="+",help="path to one or more Gadget2 snapshots to display")
	parser.add_argument("-e","--enable-mpi",dest="enable_mpi",action="store_true",default=False,help="enables the import of mpi4py (can cause problems on some systems)")

	args = parser.parse_args(args)
	
	if args.filename is None:
		parser.print_help()
		sys.exit(0)

	if not args.enable_mpi:
		sys.modules["mpi4py"]=None

	try:
		
		for filename in args.filename:
			snap = Gadget2Snapshot.open(filename)
			
			print("")
			print(filename+":")
			print("")
			for item in snap.header.keys():
				print("{0} = {1}".format(item,snap.header[item]))
			print("")
			
			snap.close()

	except Exception,e:
		print(e)
