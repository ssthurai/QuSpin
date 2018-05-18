import numpy as _np


def consolidate_static(static_list):
	eps = 10 * _np.finfo(_np.float64).eps

	static_dict={}
	for opstr,bonds in static_list:
		if opstr not in static_dict:
			static_dict[opstr] = {}

		for bond in bonds:
			J = bond[0]
			indx = tuple(bond[1:])
			if indx in static_dict[opstr]:
				static_dict[opstr][indx] += J
			else:
				static_dict[opstr][indx] = J

				

	static_list = []
	for opstr,opstr_dict in static_dict.items():
		for indx,J in opstr_dict.items():
			if _np.abs(J) > eps:
				static_list.append((opstr,indx,J))


	return static_list


def consolidate_dynamic(dynamic_list):
	eps = 10 * _np.finfo(_np.float64).eps
	
	dynamic_dict={}
	for opstr,bonds,f,f_args in dynamic_list:
		f_args = tuple(f_args)
		if (opstr,f,f_args) not in dynamic_dict:
			dynamic_dict[(opstr,f,f_args)] = {}

		for bond in bonds:
			J = bond[0]
			indx = tuple(bond[1:])
			if indx in dynamic_dict[(opstr,f,f_args)]:
				dynamic_dict[(opstr,f,f_args)][indx] += J
			else:
				dynamic_dict[(opstr,f,f_args)][indx] = J


	dynamic_list = []
	for (opstr,f,f_args),opstr_dict in dynamic_dict.items():
		for indx,J in opstr_dict.items():
			if _np.abs(J) > eps:
				dynamic_list.append((opstr,indx,J,f,f_args))


	return dynamic_list

def check_static(sub_list):
	"""Checks format of static list. """
	if (type(sub_list) in [list,tuple]):
 		if (len(sub_list) == 2):
			if type(sub_list[0]) is not str: raise TypeError('expecting string type for opstr')
			if type(sub_list[1]) in [list,tuple]:
				for sub_sub_list in sub_list[1]:
					if (type(sub_sub_list) in [list,tuple]) and (len(sub_sub_list) > 0):
						for element in sub_sub_list:
							if not _np.isscalar(element): raise TypeError('expecting scalar elements of indx')
					else: raise TypeError('expecting list for indx') 
			else: raise TypeError('expecting a list of one or more indx')
			return 0
		elif (len(sub_list) == 3):
			if _np.array(sub_list[0]).ndim>0: raise TypeError("expecting scalar for coupling")
			if type(sub_list[1]) is not str: raise TypeError('expecting string type for opstr')
			if (type(sub_sub[2]) in [list,tuple]) and (len(sub_sub_list) > 0):
				for element in sub_sub_list:
					if not _np.isscalar(element): raise TypeError('expecting scalar elements of indx')
			else: raise TypeError('expecting list for index, even for single site operators.') 
			return 1
		else:
			return 2
	else: 
		return 2
	

def check_dynamic(sub_list):
	"""Checks format of dynamic list. """
	if (type(sub_list) in [list,tuple]):
		if (len(sub_list) == 5):
			if _np.array(sub_list[0]).ndim>0: raise TypeError("expecting scalar for coupling")
			if type(sub_list[1]) is not str: raise TypeError('expecting string type for opstr')
			if (type(sub_sub[2]) in [list,tuple]) and (len(sub_sub_list) > 0):
				for element in sub_sub_list:
					if not _np.isscalar(element): raise TypeError('expecting scalar elements of indx')
			else: raise TypeError('expecting a list of one or more indx')
			if not hasattr(sub_list[3],"__call__"): raise TypeError('expecting callable object for driving function')
			if type(sub_list[4]) not in [list,tuple]: raise TypeError('expecting list for function arguments')
			return 0
		if (len(sub_list) == 4):
			if type(sub_list[0]) is not str: raise TypeError('expecting string type for opstr')
			if type(sub_list[1]) in [list,tuple]:
				for sub_sub_list in sub_list[1]:
					if (type(sub_sub_list) in [list,tuple]) and (len(sub_sub_list) > 0):
						for element in sub_sub_list:
							if not _np.isscalar(element): raise TypeError('expecting scalar elements of indx')
					else: raise TypeError('expecting list for indx') 
			else: raise TypeError('expecting a list of one or more indx')
			if not hasattr(sub_list[2],"__call__"): raise TypeError('expecting callable object for driving function')
			if type(sub_list[3]) not in [list,tuple]: raise TypeError('expecting list for function arguments')
			return 1
		elif (len(sub_list) == 3):
			if not hasattr(sub_list[1],"__call__"): raise TypeError('expecting callable object for driving function')
			if type(sub_list[2]) not in [list,tuple]: raise TypeError('expecting list for function arguments')
			return 2
		elif (len(sub_list) == 2):
			if not hasattr(sub_list[1],"__call__"): raise TypeError('expecting callable object for driving function')
			return 2
	else:
		raise TypeError('expecting list with object, driving function, and function arguments')





def _process_hamiltonian_lists(static_list,dynamic_list):
	eps = 10 * _np.finfo(_np.float64).eps

	static_opstr={}
	static_other_list=[]

	if type(static_list) in [list,tuple]:
		for ele in static_list:
			i = check_static(ele)
			if i==0:
				J,opstr,index = ele
				index = tuple(index)
				if opstr not in static_opstr:
					static_opstr[opstr] = {}

				if index in static_opstr[opstr]:
					static_opstr[opstr][index] += J
				else:
					static_opstr[opstr][index] = J

			elif i==1:
				opstr,bonds = ele
				if opstr not in static_opstr:
					static_opstr[opstr] = {}

				for bond in bonds:
					J,index = bond[0],tuple(bond[1:])
					if index in static_opstr[opstr]:
						static_opstr[opstr][index] += J
					else:
						static_opstr[opstr][index] = J
			else:
				static_other_list.append(ele)
	else: 
		raise TypeError('expecting list/tuple of lists/tuples containing opstr and list of indx')


	static_opstr_list = []
	for opstr,indices in static_opstr.items():
		for index,J in indices.items():
			if J > eps:
				static_opstr_list.append([J,opstr,index])


	dynamic_opstr={}
	dynamic_other_list=[]
	if type(dynamic_list) in [list,tuple]:
		for ele in dynamic_list:
			i = check_dynamic(ele)
			if i==0:
				J,opstr,index,f,f_args = ele
				index = tuple(index)
				key = (f,f_args,opstr)

				if key not in dynamic_opstr:
					dynamic_opstr[key] = {}

				if index in dynamic_opstr[key]:
					dynamic_opstr[key][index] += J
				else:
					dynamic_opstr[key][index] = J

			elif i==1:
				opstr,bonds,f,f_args = ele
				key = (f,f_args,opstr)

				if key not in dynamic_opstr:
					dynamic_opstr[key] = {}

				for bond in bonds:
					J,index = bond[0],tuple(bond[1:])
					if index in dynamic_opstr[key]:
						dynamic_opstr[key][index] += J
					else:
						dynamic_opstr[key][index] = J
			else: 
				dynamic_other_list.append(ele)					
	else: 
		raise TypeError('expecting list/tuple of lists/tuples containing opstr and list of indx, functions, and function args')



	dynamic_opstr_list = []
	for (f,f_args,opstr),indices in dynamic_opstr.items():
		for index,J in indices.items():
			if J > eps:
				dynamic_opstr_list.append([J,f,f_args,opstr,index])


	return static_opstr_list,static_other_list,dynamic_opstr_list,dynamic_other_list




