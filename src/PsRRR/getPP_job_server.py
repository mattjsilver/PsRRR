# PATHWAYS SPARSE REDUCED-RANK REGRESSION
#
# Distributed under GNU general public licence (see COPYRIGHT.txt) 
# Copyright (C) 2012 Matt Silver
#
# WARRANTY DISCLAIMER
# Please note there is no warranty for the program, to the extent permitted by applicable law.
# Except when otherwise stated in writing the copyright holders and/or other parties provide the program
# "as is" without warranty of any kind, either expressed or implied, including, but not limited to, 
# the implied warranties of merchantability and fitness for a particular purpose. The entire risk
# as to the quality and performance of the program is with you. Should the program
# prove defective, you assume the cost of all necessary servicing, repair or correction.
#
# Please send comments or bugs to g.montana@imperial.ac.uk
#
#**************************************************************************************************
'''
launch parallel python job server
'''

import pp, sys, time, os

class Server(pp.Server):
	'''
	create job_server instance
	'''


	def __init__(self,params):
		'''
		Constructor
		'''

		ppservers = params.ppServers # tuple of all parallel python servers to connect with

		# Create jobserver with automatically detected number of workers
 		pp.Server.__init__(self,ncpus=params.nLocalCpus, ppservers=ppservers, restart=False)

		# no response from remote servers - try launching ppserver.py
		if params.nLocalCpus: # running ppserver locally (used for testing)

			try:
				print 'launching ppserver.py locally'
				cmd = os.system(params.PsRRR_ROOT_DIR + 'src/PsRRR/ppserver.py -t 60 -p ' + \
					str(params.ppPort) + ' -w ' + str(params.workersPerServer) + ' &')
			except:
				print 'unable to launch local ppserver.py... exiting'
				sys.exit()
	
		else: # running ppserver on remote machine(s)

			try:
				print 'launching ppserver.py on remote servers'
				for server_idx in range(len(params.ppServers)):
					ppServer_name = params.ppServers[server_idx]
					# launch ppserver.py on all remote machines
					print 'launching server:', ppServer_name
					cmd = os.system('ssh ' + ppServer_name + ' ' + params.PsRRR_ROOT_DIR + \
						'src/PsRRR/ppserver.py -t 60 -p ' + str(params.ppPort) + \
						' -w ' + str(params.workersPerServer) + ' &')
					
				self.nWorkers = self.get_active_nodes() # number of available workers
			except:
				print 'unable to launch ppserver.py... exiting'
				sys.exit()

			print 'waiting for response from remote servers...'
			time.sleep(10) # pause for remote ppservers to start

		if params.nLocalCpus:
			self.nWorkers = self.get_active_nodes().values()[0] # number of available workers
		else:
			self.nWorkers = sum(self.get_active_nodes().values())

#		import interactiveShell as ipy; ipy.ipshell() # start ipython shell for debugging

		print 'active nodes for parallel processing:',	self.get_active_nodes()
		nActiveNodes = len(self.get_active_nodes())-1 # number of servers available
		print "Starting parallel python with", nActiveNodes, 'remote servers, and', self.get_ncpus(), "local workers"

#  	import interactiveShell as ipy; ipy.ipshell() # start ipython shell for debugging