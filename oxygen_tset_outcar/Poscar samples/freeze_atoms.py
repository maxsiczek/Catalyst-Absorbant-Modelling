from ase.io import vasp
from ase.constraints import FixAtoms
import os


for i in range(0,20):
    os.chdir(str(i))
    atoms_obj=vasp.read_vasp('POSCAR')
    c = FixAtoms(indices=[atom.index for atom in atoms_obj if atoms_obj.get_positions()[atom.index][2] < 11])
    atoms_obj.set_constraint(c)
    vasp.write_vasp('POSCAR',atoms_obj,sort=True)
    os.chdir('..')


def prepPoscar():
    def replace_line(file_name, line_num, text):
        lines = open(file_name, 'r').readlines()
        lines[line_num] = text
        out = open(file_name, 'w')
        out.writelines(lines)
        out.close()

    for i in range(0, 20):
        os.chdir(str(i))
        with open(r"POSCAR", 'r+') as fp:
            # read an store all lines into list
            lines = fp.readlines()
            # move file pointer to the beginning of a file
            fp.seek(0)
            # truncate the file
            fp.truncate()
            fp.writelines(lines[:-8])
            # start writing lines
            # iterate line and line number
        replace_line("POSCAR",0,'Au  O Pd\n')
        replace_line("POSCAR", 5, ' Au  O   Pd\n')
        replace_line("POSCAR", 6, '  24   8  24\n')
        os.chdir('..')

# prepPoscar()
