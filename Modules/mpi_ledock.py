import os
import sys
import mpi4py.MPI as MPI


def main(currentPath, result_path, ledockin_path):
    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    pdbname = os.path.basename(ledockin_path)
    ledockinpath = os.path.join(ledockin_path, pdbname + '_ledock' + str(comm_rank) + '.in')
    os.system('sh ' + os.path.join(currentPath, 'dock.sh ') + result_path + ledockin_path.rpartition('ledockin')[2] + ' ' + ledockinpath)


if __name__ == "__main__":
    currentPath = sys.argv[1]
    result_path = sys.argv[2]
    ledockin_path = sys.argv[3]

    main(currentPath, result_path, ledockin_path)
