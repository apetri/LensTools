from ..simulations import Gadget2Snapshot
import sys
import argparse

def main(args=None):

	parser = argparse.ArgumentParser()
	parser.add_argument("filename",nargs="+",help="path to one or more Gadget2 snapshots to display")

	args = parser.parse_args(args)
	
	if args.filename is None:
		parser.print_help()
		sys.exit(0)

	try:
		
		for filename in args.filename:
			snap = Gadget2Snapshot.open(filename)
			
			print("")
			print(filename+":")
			print("")
			keys = snap.header.keys()
			keys.sort()
			for item in keys:
				print("{0} = {1}".format(item,snap.header[item]))
			print("")
			
			snap.close()

	except Exception,e:
		print(e)
