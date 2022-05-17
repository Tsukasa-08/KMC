#!/usr/bin/env python
# coding: utf-8

# In[20]:


import pymatgen as mg
import argparse


# In[21]:


parser = argparse.ArgumentParser(description='make supercell by input factor')


# In[22]:


parser.add_argument('-f', '--unit_cell_file', type=str,                         help='Structure file name in the cif or POSCAR format                                containing only sites of diffusion species                                with symmetry information.')
parser.add_argument('-s', '--scaling_factor', type=int,                         help='A number, which simply scales all lattice vectors                                by the same factor.')

args = parser.parse_args()


# In[23]:


structure = mg.Structure.from_file(args.unit_cell_file)


# In[ ]:





# In[ ]:


structure.make_supercell(scaling_matrix=args.scaling_factor)


# In[24]:


structure.to(fmt='poscar', filename='supercell_{}.vasp'.format(args.scaling_factor))


# In[ ]:




