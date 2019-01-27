import nexusformat.nexus as nx
import pandas as pd
import matplotlib

pd.set_option('display.max_rows',8)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width',999)

scandata_field_list = ['/entry1/measurement', '/entry1/plotted']
scan_command_field_list = ['/entry1/scan_command']


class pdnx(pd.DataFrame): 
    '''
    nexusformat wrapper: tries to create dataframe from default data
    whole nexus file is under .nx attribute
    e.g. 
    n=pdnx(p % 633777)  open file for scan 633777 (p is filename/format specifier)
    n                   display pandas dataframe
    n.plot()            plot pandas dataframe
    n.plot('idgap','ic1monitor')	pandas plot with selected x and y collumns
    n.plt('idgap','ic1monitor')		same but with pdnx defaults (title etc)	
    n.nx                nexus tree
    print n.nx.tree     print nexus tree
    n.find('chi')	find 'chi' key(s) in tree and display value(s) (n.find() for all)
    n.findkeys('chi')	return list of key value lists for key 'chi'
    n.pruned_tree(n)    return nexus tree up to n levels deep
    n.nx.plot()         default nexus plot
    for i in range(633777, 633779):print pdnx(p % i).scan     print scan string for range of scans

    n['newkey'] = n.nx.entry1.before_scan.myval 	as long as 'newkey' is new then this pads out a new scan column with myval

    '''

    def __init__(self,  filestr, scandata_field_list = scandata_field_list, scan_command_field_list = scan_command_field_list):
        try:
            _nx = nx.nxload(filestr,'r')

        except:
            print("=== Error loading file %s" % filestr)
            return
        
        _load_dataframe_success = False
        for scandata_field in scandata_field_list:
            try:
                nx_scan_dict = dict(_nx[scandata_field])
                for key in nx_scan_dict.keys():
                    try:
                        nx_scan_dict[key] = nx_scan_dict[key].nxdata.flatten()
                    except:
                        pass
                pd.DataFrame.__init__(self, nx_scan_dict)
                _load_dataframe_success = True
                break
            except:
                pass

        if not _load_dataframe_success:
            #print('=== Failed to create DataFrame from data - create empty DataFrame')
            pd.DataFrame.__init__(self)

        setattr(self,'nx',_nx)
        
        for scan_command in scan_command_field_list:	
            try:
                setattr(self, 'scan', filestr+'\n'+_nx[scan_command].nxdata)
                break
            except:
                pass


    def plt(self, *args, **kwargs):
        kwargs.setdefault('title', self.scan)
        kwargs.setdefault('grid', True)
        self.plot(*args, **kwargs)


    def _list_to_dot_sep_string(self, lst):
        outstr = ''
        for item in lst:
            outstr += '.' + str(item)
        return outstr

    def _find_key(self, tree, key, previous_keys=[]):
        global _keylist
        try:
            for keyval in tree.keys():
                if keyval == key or key == '':
                    _keylist += [previous_keys + [keyval]]
                self._find_key(tree[keyval], key, previous_keys = previous_keys + [keyval])
        except:
            pass

    def findkeys(self, keystring):
        'Return list of key sequences (lists) that end with keystring'    
        global _keylist
        _keylist=[]
        self._find_key(self.nx, keystring)
        return _keylist

    def find(self, keystring=''):
        'Return nexus fields and values for keystring'
        for key_sequence in self.findkeys(keystring):
            obj = self.nx
            for key in key_sequence:
                obj = obj[key]
            print '.nx' + self._list_to_dot_sep_string(key_sequence) + ' : \t', obj #change syntax for Python3

    def pruned_tree(self, depth):
        'Print pruned tree'
        allfieldlist = self.findkeys('')
        previous = []
        for fieldlist in allfieldlist:
            fieldshort = fieldlist[:depth]
            if fieldshort != previous:
                print(self._list_to_dot_sep_string(fieldshort))
                previous = fieldshort





def vec2mat(vecx, vecy, vecz, n_inner=None):
    #matx, maty, matz = vec2mat(vecx, vecy, vecz, n_inner=None)
    #convert vectors from 2D scan to matrices
    #vecx,y,z: arrays (any dimension) or lists
    #matx,y,z: 2D arrays
    #n_inner: Number of points in inner loop - calculated if not specified
    #Arrays are truncated if the size doesn't match the required shape

    vx = np.array(vecx[:]); vy = np.array(vecy[:]); vz = np.array(vecz[:]) #get inputs in standard form
    if n_inner == None:   #calculate number in inner loop by looking for jumps
        jumps = np.abs(np.diff(vx) * np.diff(vy))
        n_inne = matplotlib.mlab.find(jumps>np.mean(jumps))[0] + 1

    n_outer = len(vx)/n_inner
    #reshape matrices
    matx = vx[0:n_inner * n_outer].reshape(n_outer,n_inner)
    maty = vy[0:n_inner * n_outer].reshape(n_outer,n_inner)
    matz = vz[0:n_inner * n_outer].reshape(n_outer,n_inner)

    return matx, maty, matz





