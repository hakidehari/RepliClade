import sys
from util.file_util import FileStream
from util.genbank_connector import GenBankConnector
from sim.simulator import Simulator


simulator = Simulator()


if __name__ == '__main__':
    
   simulator.run_simulation()