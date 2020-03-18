from __future__ import print_function, absolute_import, division
import _alegendre
import f90wrap.runtime
import logging

class Cone(f90wrap.runtime.FortranModule):
    """
    Module cone
    
    
    Defined at cone.fpp lines 5-65
    
    """
    @staticmethod
    def alegendre_array_p_eval(dnu, dmu, t_array, res):
        """
        alegendre_array_p_eval(dnu, dmu, t_array, res)
        
        
        Defined at cone.fpp lines 11-28
        
        Parameters
        ----------
        dnu : unknown
        dmu : unknown
        t_array : unknown array
        res : unknown array
        
        """
        _alegendre.f90wrap_alegendre_array_p_eval(dnu=dnu, dmu=dmu, t_array=t_array, \
            res=res)
    
    @staticmethod
    def alegendre_array_proots(dnu, dmu, roots):
        """
        alegendre_array_proots(dnu, dmu, roots)
        
        
        Defined at cone.fpp lines 30-44
        
        Parameters
        ----------
        dnu : unknown
        dmu : unknown
        roots : unknown array
        
        """
        _alegendre.f90wrap_alegendre_array_proots(dnu=dnu, dmu=dmu, roots=roots)
    
    @staticmethod
    def alegendre_array_jacobi(dmu, n, t_array, wht_array):
        """
        alegendre_array_jacobi(dmu, n, t_array, wht_array)
        
        
        Defined at cone.fpp lines 46-65
        
        Parameters
        ----------
        dmu : unknown
        n : int
        t_array : unknown array
        wht_array : unknown array
        
        """
        _alegendre.f90wrap_alegendre_array_jacobi(dmu=dmu, n=n, t_array=t_array, \
            wht_array=wht_array)
    
    _dt_array_initialisers = []
    

cone = Cone()

class Alegendreeval(f90wrap.runtime.FortranModule):
    """
    Module alegendreeval
    
    
    Defined at alegendre_eval.fpp lines 115-2834
    
    """
    @f90wrap.runtime.register_class("alegendre.alegendre_expansion_data")
    class alegendre_expansion_data(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=alegendre_expansion_data)
        
        
        Defined at alegendre_eval.fpp lines 124-180
        
        """
        def __init__(self, handle=None):
            """
            self = Alegendre_Expansion_Data()
            
            
            Defined at alegendre_eval.fpp lines 124-180
            
            
            Returns
            -------
            this : Alegendre_Expansion_Data
            	Object to be constructed
            
            
            Automatically generated constructor for alegendre_expansion_data
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _alegendre.f90wrap_alegendre_expansion_data_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Alegendre_Expansion_Data
            
            
            Defined at alegendre_eval.fpp lines 124-180
            
            Parameters
            ----------
            this : Alegendre_Expansion_Data
            	Object to be destructed
            
            
            Automatically generated destructor for alegendre_expansion_data
            """
            if self._alloc:
                _alegendre.f90wrap_alegendre_expansion_data_finalise(this=self._handle)
        
        @property
        def ifbetas(self):
            """
            Element ifbetas ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 126
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__ifbetas(self._handle)
        
        @ifbetas.setter
        def ifbetas(self, ifbetas):
            _alegendre.f90wrap_alegendre_expansion_data__set__ifbetas(self._handle, ifbetas)
        
        @property
        def ifover(self):
            """
            Element ifover ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 126
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__ifover(self._handle)
        
        @ifover.setter
        def ifover(self, ifover):
            _alegendre.f90wrap_alegendre_expansion_data__set__ifover(self._handle, ifover)
        
        @property
        def ifsmalldmu(self):
            """
            Element ifsmalldmu ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 126
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__ifsmalldmu(self._handle)
        
        @ifsmalldmu.setter
        def ifsmalldmu(self, ifsmalldmu):
            _alegendre.f90wrap_alegendre_expansion_data__set__ifsmalldmu(self._handle, \
                ifsmalldmu)
        
        @property
        def ifinverse(self):
            """
            Element ifinverse ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 126
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__ifinverse(self._handle)
        
        @ifinverse.setter
        def ifinverse(self, ifinverse):
            _alegendre.f90wrap_alegendre_expansion_data__set__ifinverse(self._handle, \
                ifinverse)
        
        @property
        def dnu1(self):
            """
            Element dnu1 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 127
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__dnu1(self._handle)
        
        @dnu1.setter
        def dnu1(self, dnu1):
            _alegendre.f90wrap_alegendre_expansion_data__set__dnu1(self._handle, dnu1)
        
        @property
        def dnu2(self):
            """
            Element dnu2 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 127
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__dnu2(self._handle)
        
        @dnu2.setter
        def dnu2(self, dnu2):
            _alegendre.f90wrap_alegendre_expansion_data__set__dnu2(self._handle, dnu2)
        
        @property
        def ncoefs(self):
            """
            Element ncoefs ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 128
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__ncoefs(self._handle)
        
        @ncoefs.setter
        def ncoefs(self, ncoefs):
            _alegendre.f90wrap_alegendre_expansion_data__set__ncoefs(self._handle, ncoefs)
        
        @property
        def k(self):
            """
            Element k ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 128
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__k(self._handle)
        
        @k.setter
        def k(self, k):
            _alegendre.f90wrap_alegendre_expansion_data__set__k(self._handle, k)
        
        @property
        def nintsab(self):
            """
            Element nintsab ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 128
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__nintsab(self._handle)
        
        @nintsab.setter
        def nintsab(self, nintsab):
            _alegendre.f90wrap_alegendre_expansion_data__set__nintsab(self._handle, nintsab)
        
        @property
        def nintscd(self):
            """
            Element nintscd ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 128
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__nintscd(self._handle)
        
        @nintscd.setter
        def nintscd(self, nintscd):
            _alegendre.f90wrap_alegendre_expansion_data__set__nintscd(self._handle, nintscd)
        
        @property
        def nintsef(self):
            """
            Element nintsef ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 128
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__nintsef(self._handle)
        
        @nintsef.setter
        def nintsef(self, nintsef):
            _alegendre.f90wrap_alegendre_expansion_data__set__nintsef(self._handle, nintsef)
        
        @property
        def ab(self):
            """
            Element ab ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 129
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__ab(self._handle)
            if array_handle in self._arrays:
                ab = self._arrays[array_handle]
            else:
                ab = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__ab)
                self._arrays[array_handle] = ab
            return ab
        
        @ab.setter
        def ab(self, ab):
            self.ab[...] = ab
        
        @property
        def cd(self):
            """
            Element cd ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 129
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__cd(self._handle)
            if array_handle in self._arrays:
                cd = self._arrays[array_handle]
            else:
                cd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__cd)
                self._arrays[array_handle] = cd
            return cd
        
        @cd.setter
        def cd(self, cd):
            self.cd[...] = cd
        
        @property
        def ef(self):
            """
            Element ef ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 129
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__ef(self._handle)
            if array_handle in self._arrays:
                ef = self._arrays[array_handle]
            else:
                ef = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__ef)
                self._arrays[array_handle] = ef
            return ef
        
        @ef.setter
        def ef(self, ef):
            self.ef[...] = ef
        
        @property
        def ncoefsalpha(self):
            """
            Element ncoefsalpha ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 135
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__ncoefsalpha(self._handle)
        
        @ncoefsalpha.setter
        def ncoefsalpha(self, ncoefsalpha):
            _alegendre.f90wrap_alegendre_expansion_data__set__ncoefsalpha(self._handle, \
                ncoefsalpha)
        
        @property
        def ncoefsalphap(self):
            """
            Element ncoefsalphap ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 135
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__ncoefsalphap(self._handle)
        
        @ncoefsalphap.setter
        def ncoefsalphap(self, ncoefsalphap):
            _alegendre.f90wrap_alegendre_expansion_data__set__ncoefsalphap(self._handle, \
                ncoefsalphap)
        
        @property
        def coefsalpha(self):
            """
            Element coefsalpha ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 136
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalpha(self._handle)
            if array_handle in self._arrays:
                coefsalpha = self._arrays[array_handle]
            else:
                coefsalpha = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalpha)
                self._arrays[array_handle] = coefsalpha
            return coefsalpha
        
        @coefsalpha.setter
        def coefsalpha(self, coefsalpha):
            self.coefsalpha[...] = coefsalpha
        
        @property
        def coefsalphap(self):
            """
            Element coefsalphap ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 136
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphap(self._handle)
            if array_handle in self._arrays:
                coefsalphap = self._arrays[array_handle]
            else:
                coefsalphap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphap)
                self._arrays[array_handle] = coefsalphap
            return coefsalphap
        
        @coefsalphap.setter
        def coefsalphap(self, coefsalphap):
            self.coefsalphap[...] = coefsalphap
        
        @property
        def iptrsalpha(self):
            """
            Element iptrsalpha ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 137
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalpha(self._handle)
            if array_handle in self._arrays:
                iptrsalpha = self._arrays[array_handle]
            else:
                iptrsalpha = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalpha)
                self._arrays[array_handle] = iptrsalpha
            return iptrsalpha
        
        @iptrsalpha.setter
        def iptrsalpha(self, iptrsalpha):
            self.iptrsalpha[...] = iptrsalpha
        
        @property
        def iptrsalphap(self):
            """
            Element iptrsalphap ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 137
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalphap(self._handle)
            if array_handle in self._arrays:
                iptrsalphap = self._arrays[array_handle]
            else:
                iptrsalphap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalphap)
                self._arrays[array_handle] = iptrsalphap
            return iptrsalphap
        
        @iptrsalphap.setter
        def iptrsalphap(self, iptrsalphap):
            self.iptrsalphap[...] = iptrsalphap
        
        @property
        def ncoefsalphapp(self):
            """
            Element ncoefsalphapp ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 143
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__ncoefsalphapp(self._handle)
        
        @ncoefsalphapp.setter
        def ncoefsalphapp(self, ncoefsalphapp):
            _alegendre.f90wrap_alegendre_expansion_data__set__ncoefsalphapp(self._handle, \
                ncoefsalphapp)
        
        @property
        def coefsalphapp(self):
            """
            Element coefsalphapp ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 144
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphapp(self._handle)
            if array_handle in self._arrays:
                coefsalphapp = self._arrays[array_handle]
            else:
                coefsalphapp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphapp)
                self._arrays[array_handle] = coefsalphapp
            return coefsalphapp
        
        @coefsalphapp.setter
        def coefsalphapp(self, coefsalphapp):
            self.coefsalphapp[...] = coefsalphapp
        
        @property
        def iptrsalphapp(self):
            """
            Element iptrsalphapp ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 145
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalphapp(self._handle)
            if array_handle in self._arrays:
                iptrsalphapp = self._arrays[array_handle]
            else:
                iptrsalphapp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalphapp)
                self._arrays[array_handle] = iptrsalphapp
            return iptrsalphapp
        
        @iptrsalphapp.setter
        def iptrsalphapp(self, iptrsalphapp):
            self.iptrsalphapp[...] = iptrsalphapp
        
        @property
        def ncoefsbeta1(self):
            """
            Element ncoefsbeta1 ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 151
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__ncoefsbeta1(self._handle)
        
        @ncoefsbeta1.setter
        def ncoefsbeta1(self, ncoefsbeta1):
            _alegendre.f90wrap_alegendre_expansion_data__set__ncoefsbeta1(self._handle, \
                ncoefsbeta1)
        
        @property
        def ncoefsbeta2(self):
            """
            Element ncoefsbeta2 ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 151
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__ncoefsbeta2(self._handle)
        
        @ncoefsbeta2.setter
        def ncoefsbeta2(self, ncoefsbeta2):
            _alegendre.f90wrap_alegendre_expansion_data__set__ncoefsbeta2(self._handle, \
                ncoefsbeta2)
        
        @property
        def nintsabb(self):
            """
            Element nintsabb ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 151
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__nintsabb(self._handle)
        
        @nintsabb.setter
        def nintsabb(self, nintsabb):
            _alegendre.f90wrap_alegendre_expansion_data__set__nintsabb(self._handle, \
                nintsabb)
        
        @property
        def coefsbeta1(self):
            """
            Element coefsbeta1 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 152
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta1(self._handle)
            if array_handle in self._arrays:
                coefsbeta1 = self._arrays[array_handle]
            else:
                coefsbeta1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta1)
                self._arrays[array_handle] = coefsbeta1
            return coefsbeta1
        
        @coefsbeta1.setter
        def coefsbeta1(self, coefsbeta1):
            self.coefsbeta1[...] = coefsbeta1
        
        @property
        def coefsbeta2(self):
            """
            Element coefsbeta2 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 152
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta2(self._handle)
            if array_handle in self._arrays:
                coefsbeta2 = self._arrays[array_handle]
            else:
                coefsbeta2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta2)
                self._arrays[array_handle] = coefsbeta2
            return coefsbeta2
        
        @coefsbeta2.setter
        def coefsbeta2(self, coefsbeta2):
            self.coefsbeta2[...] = coefsbeta2
        
        @property
        def abb(self):
            """
            Element abb ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 152
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__abb(self._handle)
            if array_handle in self._arrays:
                abb = self._arrays[array_handle]
            else:
                abb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__abb)
                self._arrays[array_handle] = abb
            return abb
        
        @abb.setter
        def abb(self, abb):
            self.abb[...] = abb
        
        @property
        def iptrsbeta1(self):
            """
            Element iptrsbeta1 ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 153
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__iptrsbeta1(self._handle)
            if array_handle in self._arrays:
                iptrsbeta1 = self._arrays[array_handle]
            else:
                iptrsbeta1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__iptrsbeta1)
                self._arrays[array_handle] = iptrsbeta1
            return iptrsbeta1
        
        @iptrsbeta1.setter
        def iptrsbeta1(self, iptrsbeta1):
            self.iptrsbeta1[...] = iptrsbeta1
        
        @property
        def iptrsbeta2(self):
            """
            Element iptrsbeta2 ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 153
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__iptrsbeta2(self._handle)
            if array_handle in self._arrays:
                iptrsbeta2 = self._arrays[array_handle]
            else:
                iptrsbeta2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__iptrsbeta2)
                self._arrays[array_handle] = iptrsbeta2
            return iptrsbeta2
        
        @iptrsbeta2.setter
        def iptrsbeta2(self, iptrsbeta2):
            self.iptrsbeta2[...] = iptrsbeta2
        
        @property
        def ncoefsalphainv(self):
            """
            Element ncoefsalphainv ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 159
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__ncoefsalphainv(self._handle)
        
        @ncoefsalphainv.setter
        def ncoefsalphainv(self, ncoefsalphainv):
            _alegendre.f90wrap_alegendre_expansion_data__set__ncoefsalphainv(self._handle, \
                ncoefsalphainv)
        
        @property
        def nintsinv(self):
            """
            Element nintsinv ftype=integer                           pytype=int
            
            
            Defined at alegendre_eval.fpp line 159
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__nintsinv(self._handle)
        
        @nintsinv.setter
        def nintsinv(self, nintsinv):
            _alegendre.f90wrap_alegendre_expansion_data__set__nintsinv(self._handle, \
                nintsinv)
        
        @property
        def coefsalphainv(self):
            """
            Element coefsalphainv ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 160
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphainv(self._handle)
            if array_handle in self._arrays:
                coefsalphainv = self._arrays[array_handle]
            else:
                coefsalphainv = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphainv)
                self._arrays[array_handle] = coefsalphainv
            return coefsalphainv
        
        @coefsalphainv.setter
        def coefsalphainv(self, coefsalphainv):
            self.coefsalphainv[...] = coefsalphainv
        
        @property
        def abinv(self):
            """
            Element abinv ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 160
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__abinv(self._handle)
            if array_handle in self._arrays:
                abinv = self._arrays[array_handle]
            else:
                abinv = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__abinv)
                self._arrays[array_handle] = abinv
            return abinv
        
        @abinv.setter
        def abinv(self, abinv):
            self.abinv[...] = abinv
        
        @property
        def iptrsalphainv(self):
            """
            Element iptrsalphainv ftype=integer pytype=int
            
            
            Defined at alegendre_eval.fpp line 161
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalphainv(self._handle)
            if array_handle in self._arrays:
                iptrsalphainv = self._arrays[array_handle]
            else:
                iptrsalphainv = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__iptrsalphainv)
                self._arrays[array_handle] = iptrsalphainv
            return iptrsalphainv
        
        @iptrsalphainv.setter
        def iptrsalphainv(self, iptrsalphainv):
            self.iptrsalphainv[...] = iptrsalphainv
        
        @property
        def epsrequired(self):
            """
            Element epsrequired ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 167
            
            """
            return \
                _alegendre.f90wrap_alegendre_expansion_data__get__epsrequired(self._handle)
        
        @epsrequired.setter
        def epsrequired(self, epsrequired):
            _alegendre.f90wrap_alegendre_expansion_data__set__epsrequired(self._handle, \
                epsrequired)
        
        @property
        def epsphase(self):
            """
            Element epsphase ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 167
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__epsphase(self._handle)
        
        @epsphase.setter
        def epsphase(self, epsphase):
            _alegendre.f90wrap_alegendre_expansion_data__set__epsphase(self._handle, \
                epsphase)
        
        @property
        def epsdisc(self):
            """
            Element epsdisc ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 167
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__epsdisc(self._handle)
        
        @epsdisc.setter
        def epsdisc(self, epsdisc):
            _alegendre.f90wrap_alegendre_expansion_data__set__epsdisc(self._handle, epsdisc)
        
        @property
        def epscomp(self):
            """
            Element epscomp ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 167
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__epscomp(self._handle)
        
        @epscomp.setter
        def epscomp(self, epscomp):
            _alegendre.f90wrap_alegendre_expansion_data__set__epscomp(self._handle, epscomp)
        
        @property
        def dmemory(self):
            """
            Element dmemory ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 168
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__dmemory(self._handle)
        
        @dmemory.setter
        def dmemory(self, dmemory):
            _alegendre.f90wrap_alegendre_expansion_data__set__dmemory(self._handle, dmemory)
        
        @property
        def dtime(self):
            """
            Element dtime ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 168
            
            """
            return _alegendre.f90wrap_alegendre_expansion_data__get__dtime(self._handle)
        
        @dtime.setter
        def dtime(self, dtime):
            _alegendre.f90wrap_alegendre_expansion_data__set__dtime(self._handle, dtime)
        
        @property
        def coefsalpha0(self):
            """
            Element coefsalpha0 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 174
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalpha0(self._handle)
            if array_handle in self._arrays:
                coefsalpha0 = self._arrays[array_handle]
            else:
                coefsalpha0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalpha0)
                self._arrays[array_handle] = coefsalpha0
            return coefsalpha0
        
        @coefsalpha0.setter
        def coefsalpha0(self, coefsalpha0):
            self.coefsalpha0[...] = coefsalpha0
        
        @property
        def coefsalphap0(self):
            """
            Element coefsalphap0 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 174
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphap0(self._handle)
            if array_handle in self._arrays:
                coefsalphap0 = self._arrays[array_handle]
            else:
                coefsalphap0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphap0)
                self._arrays[array_handle] = coefsalphap0
            return coefsalphap0
        
        @coefsalphap0.setter
        def coefsalphap0(self, coefsalphap0):
            self.coefsalphap0[...] = coefsalphap0
        
        @property
        def coefsalphapp0(self):
            """
            Element coefsalphapp0 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 174
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphapp0(self._handle)
            if array_handle in self._arrays:
                coefsalphapp0 = self._arrays[array_handle]
            else:
                coefsalphapp0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphapp0)
                self._arrays[array_handle] = coefsalphapp0
            return coefsalphapp0
        
        @coefsalphapp0.setter
        def coefsalphapp0(self, coefsalphapp0):
            self.coefsalphapp0[...] = coefsalphapp0
        
        @property
        def coefsbeta10(self):
            """
            Element coefsbeta10 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 175
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta10(self._handle)
            if array_handle in self._arrays:
                coefsbeta10 = self._arrays[array_handle]
            else:
                coefsbeta10 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta10)
                self._arrays[array_handle] = coefsbeta10
            return coefsbeta10
        
        @coefsbeta10.setter
        def coefsbeta10(self, coefsbeta10):
            self.coefsbeta10[...] = coefsbeta10
        
        @property
        def coefsbeta20(self):
            """
            Element coefsbeta20 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 175
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta20(self._handle)
            if array_handle in self._arrays:
                coefsbeta20 = self._arrays[array_handle]
            else:
                coefsbeta20 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsbeta20)
                self._arrays[array_handle] = coefsbeta20
            return coefsbeta20
        
        @coefsbeta20.setter
        def coefsbeta20(self, coefsbeta20):
            self.coefsbeta20[...] = coefsbeta20
        
        @property
        def coefsalphainv0(self):
            """
            Element coefsalphainv0 ftype=double precision pytype=unknown
            
            
            Defined at alegendre_eval.fpp line 176
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphainv0(self._handle)
            if array_handle in self._arrays:
                coefsalphainv0 = self._arrays[array_handle]
            else:
                coefsalphainv0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_alegendre_expansion_data__array__coefsalphainv0)
                self._arrays[array_handle] = coefsalphainv0
            return coefsalphainv0
        
        @coefsalphainv0.setter
        def coefsalphainv0(self, coefsalphainv0):
            self.coefsalphainv0[...] = coefsalphainv0
        
        def __str__(self):
            ret = ['<alegendre_expansion_data>{\n']
            ret.append('    ifbetas : ')
            ret.append(repr(self.ifbetas))
            ret.append(',\n    ifover : ')
            ret.append(repr(self.ifover))
            ret.append(',\n    ifsmalldmu : ')
            ret.append(repr(self.ifsmalldmu))
            ret.append(',\n    ifinverse : ')
            ret.append(repr(self.ifinverse))
            ret.append(',\n    dnu1 : ')
            ret.append(repr(self.dnu1))
            ret.append(',\n    dnu2 : ')
            ret.append(repr(self.dnu2))
            ret.append(',\n    ncoefs : ')
            ret.append(repr(self.ncoefs))
            ret.append(',\n    k : ')
            ret.append(repr(self.k))
            ret.append(',\n    nintsab : ')
            ret.append(repr(self.nintsab))
            ret.append(',\n    nintscd : ')
            ret.append(repr(self.nintscd))
            ret.append(',\n    nintsef : ')
            ret.append(repr(self.nintsef))
            ret.append(',\n    ab : ')
            ret.append(repr(self.ab))
            ret.append(',\n    cd : ')
            ret.append(repr(self.cd))
            ret.append(',\n    ef : ')
            ret.append(repr(self.ef))
            ret.append(',\n    ncoefsalpha : ')
            ret.append(repr(self.ncoefsalpha))
            ret.append(',\n    ncoefsalphap : ')
            ret.append(repr(self.ncoefsalphap))
            ret.append(',\n    coefsalpha : ')
            ret.append(repr(self.coefsalpha))
            ret.append(',\n    coefsalphap : ')
            ret.append(repr(self.coefsalphap))
            ret.append(',\n    iptrsalpha : ')
            ret.append(repr(self.iptrsalpha))
            ret.append(',\n    iptrsalphap : ')
            ret.append(repr(self.iptrsalphap))
            ret.append(',\n    ncoefsalphapp : ')
            ret.append(repr(self.ncoefsalphapp))
            ret.append(',\n    coefsalphapp : ')
            ret.append(repr(self.coefsalphapp))
            ret.append(',\n    iptrsalphapp : ')
            ret.append(repr(self.iptrsalphapp))
            ret.append(',\n    ncoefsbeta1 : ')
            ret.append(repr(self.ncoefsbeta1))
            ret.append(',\n    ncoefsbeta2 : ')
            ret.append(repr(self.ncoefsbeta2))
            ret.append(',\n    nintsabb : ')
            ret.append(repr(self.nintsabb))
            ret.append(',\n    coefsbeta1 : ')
            ret.append(repr(self.coefsbeta1))
            ret.append(',\n    coefsbeta2 : ')
            ret.append(repr(self.coefsbeta2))
            ret.append(',\n    abb : ')
            ret.append(repr(self.abb))
            ret.append(',\n    iptrsbeta1 : ')
            ret.append(repr(self.iptrsbeta1))
            ret.append(',\n    iptrsbeta2 : ')
            ret.append(repr(self.iptrsbeta2))
            ret.append(',\n    ncoefsalphainv : ')
            ret.append(repr(self.ncoefsalphainv))
            ret.append(',\n    nintsinv : ')
            ret.append(repr(self.nintsinv))
            ret.append(',\n    coefsalphainv : ')
            ret.append(repr(self.coefsalphainv))
            ret.append(',\n    abinv : ')
            ret.append(repr(self.abinv))
            ret.append(',\n    iptrsalphainv : ')
            ret.append(repr(self.iptrsalphainv))
            ret.append(',\n    epsrequired : ')
            ret.append(repr(self.epsrequired))
            ret.append(',\n    epsphase : ')
            ret.append(repr(self.epsphase))
            ret.append(',\n    epsdisc : ')
            ret.append(repr(self.epsdisc))
            ret.append(',\n    epscomp : ')
            ret.append(repr(self.epscomp))
            ret.append(',\n    dmemory : ')
            ret.append(repr(self.dmemory))
            ret.append(',\n    dtime : ')
            ret.append(repr(self.dtime))
            ret.append(',\n    coefsalpha0 : ')
            ret.append(repr(self.coefsalpha0))
            ret.append(',\n    coefsalphap0 : ')
            ret.append(repr(self.coefsalphap0))
            ret.append(',\n    coefsalphapp0 : ')
            ret.append(repr(self.coefsalphapp0))
            ret.append(',\n    coefsbeta10 : ')
            ret.append(repr(self.coefsbeta10))
            ret.append(',\n    coefsbeta20 : ')
            ret.append(repr(self.coefsbeta20))
            ret.append(',\n    coefsalphainv0 : ')
            ret.append(repr(self.coefsalphainv0))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def alegendre_eval(dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq):
        """
        alegendre_eval(dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq)
        
        
        Defined at alegendre_eval.fpp lines 238-439
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        alpha : float
        alphader : float
        vallogp : float
        vallogq : float
        valp : float
        valq : float
        
        """
        _alegendre.f90wrap_alegendre_eval(dnu=dnu, dmu=dmu, t=t, alpha=alpha, \
            alphader=alphader, vallogp=vallogp, vallogq=vallogq, valp=valp, valq=valq)
    
    @staticmethod
    def alegendre_eval00(dnu, dmu, t, alphader):
        """
        alegendre_eval00(dnu, dmu, t, alphader)
        
        
        Defined at alegendre_eval.fpp lines 441-506
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        alphader : float
        
        """
        _alegendre.f90wrap_alegendre_eval00(dnu=dnu, dmu=dmu, t=t, alphader=alphader)
    
    @staticmethod
    def alegendre_nroots(dnu, dmu, nproots, nqroots):
        """
        alegendre_nroots(dnu, dmu, nproots, nqroots)
        
        
        Defined at alegendre_eval.fpp lines 515-560
        
        Parameters
        ----------
        dnu : float
        dmu : float
        nproots : int
        nqroots : int
        
        """
        _alegendre.f90wrap_alegendre_nroots(dnu=dnu, dmu=dmu, nproots=nproots, \
            nqroots=nqroots)
    
    @staticmethod
    def alegendre_proot(dnu, dmu, j, t):
        """
        alegendre_proot(dnu, dmu, j, t)
        
        
        Defined at alegendre_eval.fpp lines 562-610
        
        Parameters
        ----------
        dnu : float
        dmu : float
        j : int
        t : float
        
        """
        _alegendre.f90wrap_alegendre_proot(dnu=dnu, dmu=dmu, j=j, t=t)
    
    @staticmethod
    def alegendre_qroot(dnu, dmu, j, t):
        """
        alegendre_qroot(dnu, dmu, j, t)
        
        
        Defined at alegendre_eval.fpp lines 612-660
        
        Parameters
        ----------
        dnu : float
        dmu : float
        j : int
        t : float
        
        """
        _alegendre.f90wrap_alegendre_qroot(dnu=dnu, dmu=dmu, j=j, t=t)
    
    @staticmethod
    def alegendre_jacobi(dmu, n, j, t, wht):
        """
        alegendre_jacobi(dmu, n, j, t, wht)
        
        
        Defined at alegendre_eval.fpp lines 662-746
        
        Parameters
        ----------
        dmu : float
        n : int
        j : int
        t : float
        wht : float
        
        """
        _alegendre.f90wrap_alegendre_jacobi(dmu=dmu, n=n, j=j, t=t, wht=wht)
    
    @staticmethod
    def alegendre_inverse(dnu, dmu, x, t):
        """
        alegendre_inverse(dnu, dmu, x, t)
        
        
        Defined at alegendre_eval.fpp lines 748-782
        
        Parameters
        ----------
        dnu : float
        dmu : float
        x : float
        t : float
        
        """
        _alegendre.f90wrap_alegendre_inverse(dnu=dnu, dmu=dmu, x=x, t=t)
    
    @staticmethod
    def alegendre_tp(dnu, dmu, a, alpha, alphader, alphader2):
        """
        alegendre_tp(dnu, dmu, a, alpha, alphader, alphader2)
        
        
        Defined at alegendre_eval.fpp lines 791-833
        
        Parameters
        ----------
        dnu : float
        dmu : float
        a : float
        alpha : float
        alphader : float
        alphader2 : float
        
        """
        _alegendre.f90wrap_alegendre_tp(dnu=dnu, dmu=dmu, a=a, alpha=alpha, \
            alphader=alphader, alphader2=alphader2)
    
    @staticmethod
    def alegendre_eval_init(dsize):
        """
        alegendre_eval_init(dsize)
        
        
        Defined at alegendre_eval.fpp lines 842-924
        
        Parameters
        ----------
        dsize : float
        
        """
        _alegendre.f90wrap_alegendre_eval_init(dsize=dsize)
    
    @staticmethod
    def alegendre_read_expansion(iw):
        """
        expdata = alegendre_read_expansion(iw)
        
        
        Defined at alegendre_eval.fpp lines 926-1020
        
        Parameters
        ----------
        iw : int
        
        Returns
        -------
        expdata : Alegendre_Expansion_Data
        
        """
        expdata = _alegendre.f90wrap_alegendre_read_expansion(iw=iw)
        expdata = \
            f90wrap.runtime.lookup_class("alegendre.alegendre_expansion_data").from_handle(expdata)
        return expdata
    
    @staticmethod
    def alegendre_read_double_array_binary(iw, n, data):
        """
        alegendre_read_double_array_binary(iw, n, data)
        
        
        Defined at alegendre_eval.fpp lines 1022-1029
        
        Parameters
        ----------
        iw : int
        n : int
        data : unknown array
        
        """
        _alegendre.f90wrap_alegendre_read_double_array_binary(iw=iw, n=n, data=data)
    
    @staticmethod
    def alegendre_read_double_binary(iw, data):
        """
        alegendre_read_double_binary(iw, data)
        
        
        Defined at alegendre_eval.fpp lines 1031-1037
        
        Parameters
        ----------
        iw : int
        data : float
        
        """
        _alegendre.f90wrap_alegendre_read_double_binary(iw=iw, data=data)
    
    @staticmethod
    def alegendre_read_double_binary16(iw, data):
        """
        alegendre_read_double_binary16(iw, data)
        
        
        Defined at alegendre_eval.fpp lines 1039-1045
        
        Parameters
        ----------
        iw : int
        data : float
        
        """
        _alegendre.f90wrap_alegendre_read_double_binary16(iw=iw, data=data)
    
    @staticmethod
    def alegendre_read_double_array_binary16(iw, n, data):
        """
        alegendre_read_double_array_binary16(iw, n, data)
        
        
        Defined at alegendre_eval.fpp lines 1047-1063
        
        Parameters
        ----------
        iw : int
        n : int
        data : unknown array
        
        """
        _alegendre.f90wrap_alegendre_read_double_array_binary16(iw=iw, n=n, data=data)
    
    @staticmethod
    def alegendre_read_integer_array_binary(iw, n, idata):
        """
        alegendre_read_integer_array_binary(iw, n, idata)
        
        
        Defined at alegendre_eval.fpp lines 1065-1074
        
        Parameters
        ----------
        iw : int
        n : int
        idata : int array
        
        """
        _alegendre.f90wrap_alegendre_read_integer_array_binary(iw=iw, n=n, idata=idata)
    
    @staticmethod
    def alegendre_read_integer_binary(iw, idata):
        """
        alegendre_read_integer_binary(iw, idata)
        
        
        Defined at alegendre_eval.fpp lines 1076-1083
        
        Parameters
        ----------
        iw : int
        idata : int
        
        """
        _alegendre.f90wrap_alegendre_read_integer_binary(iw=iw, idata=idata)
    
    @staticmethod
    def alegendre_macdonald(dnu, dmu, t, vallogp, vallogq, valp, valq):
        """
        alegendre_macdonald(dnu, dmu, t, vallogp, vallogq, valp, valq)
        
        
        Defined at alegendre_eval.fpp lines 1092-1248
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        vallogp : float
        vallogq : float
        valp : float
        valq : float
        
        """
        _alegendre.f90wrap_alegendre_macdonald(dnu=dnu, dmu=dmu, t=t, vallogp=vallogp, \
            vallogq=vallogq, valp=valp, valq=valq)
    
    @staticmethod
    def alegendre_macdonald_logp0(dnu, dmu, t, vallogp):
        """
        alegendre_macdonald_logp0(dnu, dmu, t, vallogp)
        
        
        Defined at alegendre_eval.fpp lines 1250-1342
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        vallogp : float
        
        """
        _alegendre.f90wrap_alegendre_macdonald_logp0(dnu=dnu, dmu=dmu, t=t, \
            vallogp=vallogp)
    
    @staticmethod
    def alegendre_macdonald_logq0(dnu, dmu, t, vallogq, dsignq):
        """
        alegendre_macdonald_logq0(dnu, dmu, t, vallogq, dsignq)
        
        
        Defined at alegendre_eval.fpp lines 1344-1434
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        vallogq : float
        dsignq : float
        
        """
        _alegendre.f90wrap_alegendre_macdonald_logq0(dnu=dnu, dmu=dmu, t=t, \
            vallogq=vallogq, dsignq=dsignq)
    
    @staticmethod
    def alegendre_gamma_ratio2(dnu, dmu, val):
        """
        alegendre_gamma_ratio2(dnu, dmu, val)
        
        
        Defined at alegendre_eval.fpp lines 1436-1483
        
        Parameters
        ----------
        dnu : float
        dmu : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_gamma_ratio2(dnu=dnu, dmu=dmu, val=val)
    
    @staticmethod
    def alegendre_log_besselj(dmu, t, vallogj, dsignj, valj):
        """
        alegendre_log_besselj(dmu, t, vallogj, dsignj, valj)
        
        
        Defined at alegendre_eval.fpp lines 1485-1526
        
        Parameters
        ----------
        dmu : float
        t : float
        vallogj : float
        dsignj : float
        valj : float
        
        """
        _alegendre.f90wrap_alegendre_log_besselj(dmu=dmu, t=t, vallogj=vallogj, \
            dsignj=dsignj, valj=valj)
    
    @staticmethod
    def alegendre_log_bessely(dmu, t, vallogy, dsigny, valy):
        """
        alegendre_log_bessely(dmu, t, vallogy, dsigny, valy)
        
        
        Defined at alegendre_eval.fpp lines 1528-1576
        
        Parameters
        ----------
        dmu : float
        t : float
        vallogy : float
        dsigny : float
        valy : float
        
        """
        _alegendre.f90wrap_alegendre_log_bessely(dmu=dmu, t=t, vallogy=vallogy, \
            dsigny=dsigny, valy=valy)
    
    @staticmethod
    def alegendre_taylor(dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq, \
        ifoscillatory):
        """
        alegendre_taylor(dnu, dmu, t, alpha, alphader, vallogp, vallogq, valp, valq, \
            ifoscillatory)
        
        
        Defined at alegendre_eval.fpp lines 1586-1626
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        alpha : float
        alphader : float
        vallogp : float
        vallogq : float
        valp : float
        valq : float
        ifoscillatory : int
        
        """
        _alegendre.f90wrap_alegendre_taylor(dnu=dnu, dmu=dmu, t=t, alpha=alpha, \
            alphader=alphader, vallogp=vallogp, vallogq=vallogq, valp=valp, valq=valq, \
            ifoscillatory=ifoscillatory)
    
    @staticmethod
    def alegendre_ptaylor(dnu, dmu, t, val):
        """
        alegendre_ptaylor(dnu, dmu, t, val)
        
        
        Defined at alegendre_eval.fpp lines 1628-1688
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_ptaylor(dnu=dnu, dmu=dmu, t=t, val=val)
    
    @staticmethod
    def alegendre_qtaylor(dnu, dmu, t, val):
        """
        alegendre_qtaylor(dnu, dmu, t, val)
        
        
        Defined at alegendre_eval.fpp lines 1690-1828
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_qtaylor(dnu=dnu, dmu=dmu, t=t, val=val)
    
    @staticmethod
    def alegendre_qtaylor0(dnu, dmu, t, val):
        """
        alegendre_qtaylor0(dnu, dmu, t, val)
        
        
        Defined at alegendre_eval.fpp lines 1830-1848
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_qtaylor0(dnu=dnu, dmu=dmu, t=t, val=val)
    
    @staticmethod
    def alegendre_ptaylor_log(dnu, dmu, t, vallog, dsign):
        """
        alegendre_ptaylor_log(dnu, dmu, t, vallog, dsign)
        
        
        Defined at alegendre_eval.fpp lines 1850-1926
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        vallog : float
        dsign : float
        
        """
        _alegendre.f90wrap_alegendre_ptaylor_log(dnu=dnu, dmu=dmu, t=t, vallog=vallog, \
            dsign=dsign)
    
    @staticmethod
    def alegendre_qtaylor_log(dnu, dmu, t, vallog):
        """
        alegendre_qtaylor_log(dnu, dmu, t, vallog)
        
        
        Defined at alegendre_eval.fpp lines 1928-2069
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        vallog : float
        
        """
        _alegendre.f90wrap_alegendre_qtaylor_log(dnu=dnu, dmu=dmu, t=t, vallog=vallog)
    
    @staticmethod
    def alegendre_qtaylor0_log(dnu, dmu, t, vallog, dsign):
        """
        alegendre_qtaylor0_log(dnu, dmu, t, vallog, dsign)
        
        
        Defined at alegendre_eval.fpp lines 2071-2097
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        vallog : float
        dsign : float
        
        """
        _alegendre.f90wrap_alegendre_qtaylor0_log(dnu=dnu, dmu=dmu, t=t, vallog=vallog, \
            dsign=dsign)
    
    @staticmethod
    def alegendre_gamma(x):
        """
        alegendre_gamma = alegendre_gamma(x)
        
        
        Defined at alegendre_eval.fpp lines 2099-2121
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        alegendre_gamma : unknown
        
        """
        alegendre_gamma = _alegendre.f90wrap_alegendre_gamma(x=x)
        return alegendre_gamma
    
    @staticmethod
    def alegendre_nearright(dnu, dmu, t, aval, apval, valp, valq):
        """
        alegendre_nearright(dnu, dmu, t, aval, apval, valp, valq)
        
        
        Defined at alegendre_eval.fpp lines 2130-2220
        
        Parameters
        ----------
        dnu : float
        dmu : float
        t : float
        aval : float
        apval : float
        valp : float
        valq : float
        
        """
        _alegendre.f90wrap_alegendre_nearright(dnu=dnu, dmu=dmu, t=t, aval=aval, \
            apval=apval, valp=valp, valq=valq)
    
    @staticmethod
    def alegendre_alphap0(dnu, dmu, val):
        """
        alegendre_alphap0(dnu, dmu, val)
        
        
        Defined at alegendre_eval.fpp lines 2222-2263
        
        Parameters
        ----------
        dnu : float
        dmu : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_alphap0(dnu=dnu, dmu=dmu, val=val)
    
    @staticmethod
    def alegendre_gammratio1(x, val):
        """
        alegendre_gammratio1(x, val)
        
        
        Defined at alegendre_eval.fpp lines 2265-2311
        
        Parameters
        ----------
        x : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_gammratio1(x=x, val=val)
    
    @staticmethod
    def alegendre_expeval(self, dnu, dmu, t, valp, valq, a, b, c):
        """
        aval, apval, bval1, bval2 = alegendre_expeval(self, dnu, dmu, t, valp, valq, a, \
            b, c)
        
        
        Defined at alegendre_eval.fpp lines 2320-2408
        
        Parameters
        ----------
        expdata : Alegendre_Expansion_Data
        dnu : unknown
        dmu : float
        t : unknown
        valp : float
        valq : float
        a : float
        b : float
        c : float
        
        Returns
        -------
        aval : unknown
        apval : unknown
        bval1 : unknown
        bval2 : unknown
        
        """
        aval, apval, bval1, bval2 = \
            _alegendre.f90wrap_alegendre_expeval(expdata=self._handle, dnu=dnu, dmu=dmu, \
            t=t, valp=valp, valq=valq, a=a, b=b, c=c)
        return aval, apval, bval1, bval2
    
    @staticmethod
    def alegendre_expeval000(self, dnu, dmu, t, apval):
        """
        alegendre_expeval000(self, dnu, dmu, t, apval)
        
        
        Defined at alegendre_eval.fpp lines 2410-2450
        
        Parameters
        ----------
        expdata : Alegendre_Expansion_Data
        dnu : float
        dmu : float
        t : float
        apval : float
        
        """
        _alegendre.f90wrap_alegendre_expeval000(expdata=self._handle, dnu=dnu, dmu=dmu, \
            t=t, apval=apval)
    
    @staticmethod
    def alegendre_expeval_inverse(self, dnu, dmu, t):
        """
        ainv = alegendre_expeval_inverse(self, dnu, dmu, t)
        
        
        Defined at alegendre_eval.fpp lines 2452-2511
        
        Parameters
        ----------
        expdata : Alegendre_Expansion_Data
        dnu : unknown
        dmu : unknown
        t : unknown
        
        Returns
        -------
        ainv : unknown
        
        """
        ainv = _alegendre.f90wrap_alegendre_expeval_inverse(expdata=self._handle, \
            dnu=dnu, dmu=dmu, t=t)
        return ainv
    
    @staticmethod
    def alegendre_expeval00(self, dnu, dmu, aval, apval, appval):
        """
        alegendre_expeval00(self, dnu, dmu, aval, apval, appval)
        
        
        Defined at alegendre_eval.fpp lines 2513-2570
        
        Parameters
        ----------
        expdata : Alegendre_Expansion_Data
        dnu : float
        dmu : float
        aval : float
        apval : float
        appval : float
        
        """
        _alegendre.f90wrap_alegendre_expeval00(expdata=self._handle, dnu=dnu, dmu=dmu, \
            aval=aval, apval=apval, appval=appval)
    
    @staticmethod
    def alegendre_evalabc(ifsmalldmu, dnu, dmu, a, b, c):
        """
        alegendre_evalabc(ifsmalldmu, dnu, dmu, a, b, c)
        
        
        Defined at alegendre_eval.fpp lines 2572-2588
        
        Parameters
        ----------
        ifsmalldmu : int
        dnu : float
        dmu : float
        a : float
        b : float
        c : float
        
        """
        _alegendre.f90wrap_alegendre_evalabc(ifsmalldmu=ifsmalldmu, dnu=dnu, dmu=dmu, \
            a=a, b=b, c=c)
    
    @staticmethod
    def compute_dmu(ifsmalldmu, dnu, dmu, dmu0):
        """
        compute_dmu(ifsmalldmu, dnu, dmu, dmu0)
        
        
        Defined at alegendre_eval.fpp lines 2590-2602
        
        Parameters
        ----------
        ifsmalldmu : int
        dnu : float
        dmu : float
        dmu0 : float
        
        """
        _alegendre.f90wrap_compute_dmu(ifsmalldmu=ifsmalldmu, dnu=dnu, dmu=dmu, \
            dmu0=dmu0)
    
    @staticmethod
    def compute_dmu0(ifsmalldmu, dnu, dmu, dmu0):
        """
        compute_dmu0(ifsmalldmu, dnu, dmu, dmu0)
        
        
        Defined at alegendre_eval.fpp lines 2604-2617
        
        Parameters
        ----------
        ifsmalldmu : int
        dnu : float
        dmu : float
        dmu0 : float
        
        """
        _alegendre.f90wrap_compute_dmu0(ifsmalldmu=ifsmalldmu, dnu=dnu, dmu=dmu, \
            dmu0=dmu0)
    
    @staticmethod
    def alegendre_tensor_eval(ncoefs, coefs, iptr, a, b, c, d, e, f, x, y, z, val):
        """
        alegendre_tensor_eval(ncoefs, coefs, iptr, a, b, c, d, e, f, x, y, z, val)
        
        
        Defined at alegendre_eval.fpp lines 2619-2657
        
        Parameters
        ----------
        ncoefs : int
        coefs : unknown array
        iptr : int
        a : float
        b : float
        c : float
        d : float
        e : float
        f : float
        x : float
        y : float
        z : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_tensor_eval(ncoefs=ncoefs, coefs=coefs, iptr=iptr, \
            a=a, b=b, c=c, d=d, e=e, f=f, x=x, y=y, z=z, val=val)
    
    @staticmethod
    def alegendre_tensor_eval2(ncoefs1, coefs1, ncoefs2, coefs2, iptr1, iptr2, a, b, \
        c, d, e, f, x, y, z, val1, val2):
        """
        alegendre_tensor_eval2(ncoefs1, coefs1, ncoefs2, coefs2, iptr1, iptr2, a, b, c, \
            d, e, f, x, y, z, val1, val2)
        
        
        Defined at alegendre_eval.fpp lines 2659-2719
        
        Parameters
        ----------
        ncoefs1 : int
        coefs1 : unknown array
        ncoefs2 : int
        coefs2 : unknown array
        iptr1 : int
        iptr2 : int
        a : float
        b : float
        c : float
        d : float
        e : float
        f : float
        x : float
        y : float
        z : float
        val1 : float
        val2 : float
        
        """
        _alegendre.f90wrap_alegendre_tensor_eval2(ncoefs1=ncoefs1, coefs1=coefs1, \
            ncoefs2=ncoefs2, coefs2=coefs2, iptr1=iptr1, iptr2=iptr2, a=a, b=b, c=c, \
            d=d, e=e, f=f, x=x, y=y, z=z, val1=val1, val2=val2)
    
    @staticmethod
    def alegendre_tensor_eval2d(ncoefs, coefs, iptr0, a, b, c, d, x, y, val):
        """
        alegendre_tensor_eval2d(ncoefs, coefs, iptr0, a, b, c, d, x, y, val)
        
        
        Defined at alegendre_eval.fpp lines 2721-2752
        
        Parameters
        ----------
        ncoefs : int
        coefs : unknown array
        iptr0 : int
        a : float
        b : float
        c : float
        d : float
        x : float
        y : float
        val : float
        
        """
        _alegendre.f90wrap_alegendre_tensor_eval2d(ncoefs=ncoefs, coefs=coefs, \
            iptr0=iptr0, a=a, b=b, c=c, d=d, x=x, y=y, val=val)
    
    @staticmethod
    def alegendre_findint(nints, ab, x, int_bn, a, b):
        """
        alegendre_findint(nints, ab, x, int_bn, a, b)
        
        
        Defined at alegendre_eval.fpp lines 2754-2778
        
        Parameters
        ----------
        nints : int
        ab : unknown array
        x : unknown
        int_bn : int
        a : unknown
        b : unknown
        
        """
        _alegendre.f90wrap_alegendre_findint(nints=nints, ab=ab, x=x, int_bn=int_bn, \
            a=a, b=b)
    
    @staticmethod
    def alegendre_chebs(x, n, pols):
        """
        alegendre_chebs(x, n, pols)
        
        
        Defined at alegendre_eval.fpp lines 2780-2834
        
        Parameters
        ----------
        x : unknown
        n : int
        pols : unknown array
        
        """
        _alegendre.f90wrap_alegendre_chebs(x=x, n=n, pols=pols)
    
    _dt_array_initialisers = []
    

alegendreeval = Alegendreeval()

class Besseleval(f90wrap.runtime.FortranModule):
    """
    Module besseleval
    
    
    Defined at bessel_eval.fpp lines 86-1904
    
    """
    @f90wrap.runtime.register_class("alegendre.bessel_expansion_data")
    class bessel_expansion_data(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=bessel_expansion_data)
        
        
        Defined at bessel_eval.fpp lines 90-107
        
        """
        def __init__(self, handle=None):
            """
            self = Bessel_Expansion_Data()
            
            
            Defined at bessel_eval.fpp lines 90-107
            
            
            Returns
            -------
            this : Bessel_Expansion_Data
            	Object to be constructed
            
            
            Automatically generated constructor for bessel_expansion_data
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _alegendre.f90wrap_bessel_expansion_data_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Bessel_Expansion_Data
            
            
            Defined at bessel_eval.fpp lines 90-107
            
            Parameters
            ----------
            this : Bessel_Expansion_Data
            	Object to be destructed
            
            
            Automatically generated destructor for bessel_expansion_data
            """
            if self._alloc:
                _alegendre.f90wrap_bessel_expansion_data_finalise(this=self._handle)
        
        @property
        def ifbetas(self):
            """
            Element ifbetas ftype=integer                           pytype=int
            
            
            Defined at bessel_eval.fpp line 91
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__ifbetas(self._handle)
        
        @ifbetas.setter
        def ifbetas(self, ifbetas):
            _alegendre.f90wrap_bessel_expansion_data__set__ifbetas(self._handle, ifbetas)
        
        @property
        def ifover(self):
            """
            Element ifover ftype=integer                           pytype=int
            
            
            Defined at bessel_eval.fpp line 91
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__ifover(self._handle)
        
        @ifover.setter
        def ifover(self, ifover):
            _alegendre.f90wrap_bessel_expansion_data__set__ifover(self._handle, ifover)
        
        @property
        def ifsmall(self):
            """
            Element ifsmall ftype=integer                           pytype=int
            
            
            Defined at bessel_eval.fpp line 91
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__ifsmall(self._handle)
        
        @ifsmall.setter
        def ifsmall(self, ifsmall):
            _alegendre.f90wrap_bessel_expansion_data__set__ifsmall(self._handle, ifsmall)
        
        @property
        def dnu1(self):
            """
            Element dnu1 ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 92
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__dnu1(self._handle)
        
        @dnu1.setter
        def dnu1(self, dnu1):
            _alegendre.f90wrap_bessel_expansion_data__set__dnu1(self._handle, dnu1)
        
        @property
        def dnu2(self):
            """
            Element dnu2 ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 92
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__dnu2(self._handle)
        
        @dnu2.setter
        def dnu2(self, dnu2):
            _alegendre.f90wrap_bessel_expansion_data__set__dnu2(self._handle, dnu2)
        
        @property
        def epsrequired(self):
            """
            Element epsrequired ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 92
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__epsrequired(self._handle)
        
        @epsrequired.setter
        def epsrequired(self, epsrequired):
            _alegendre.f90wrap_bessel_expansion_data__set__epsrequired(self._handle, \
                epsrequired)
        
        @property
        def k(self):
            """
            Element k ftype=integer                           pytype=int
            
            
            Defined at bessel_eval.fpp line 93
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__k(self._handle)
        
        @k.setter
        def k(self, k):
            _alegendre.f90wrap_bessel_expansion_data__set__k(self._handle, k)
        
        @property
        def nintsab(self):
            """
            Element nintsab ftype=integer                           pytype=int
            
            
            Defined at bessel_eval.fpp line 93
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__nintsab(self._handle)
        
        @nintsab.setter
        def nintsab(self, nintsab):
            _alegendre.f90wrap_bessel_expansion_data__set__nintsab(self._handle, nintsab)
        
        @property
        def nintscd(self):
            """
            Element nintscd ftype=integer                           pytype=int
            
            
            Defined at bessel_eval.fpp line 93
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__nintscd(self._handle)
        
        @nintscd.setter
        def nintscd(self, nintscd):
            _alegendre.f90wrap_bessel_expansion_data__set__nintscd(self._handle, nintscd)
        
        @property
        def nintsef(self):
            """
            Element nintsef ftype=integer                           pytype=int
            
            
            Defined at bessel_eval.fpp line 93
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__nintsef(self._handle)
        
        @nintsef.setter
        def nintsef(self, nintsef):
            _alegendre.f90wrap_bessel_expansion_data__set__nintsef(self._handle, nintsef)
        
        @property
        def ab(self):
            """
            Element ab ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 94
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__ab(self._handle)
            if array_handle in self._arrays:
                ab = self._arrays[array_handle]
            else:
                ab = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__ab)
                self._arrays[array_handle] = ab
            return ab
        
        @ab.setter
        def ab(self, ab):
            self.ab[...] = ab
        
        @property
        def cd(self):
            """
            Element cd ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 94
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__cd(self._handle)
            if array_handle in self._arrays:
                cd = self._arrays[array_handle]
            else:
                cd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__cd)
                self._arrays[array_handle] = cd
            return cd
        
        @cd.setter
        def cd(self, cd):
            self.cd[...] = cd
        
        @property
        def ef(self):
            """
            Element ef ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 94
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__ef(self._handle)
            if array_handle in self._arrays:
                ef = self._arrays[array_handle]
            else:
                ef = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__ef)
                self._arrays[array_handle] = ef
            return ef
        
        @ef.setter
        def ef(self, ef):
            self.ef[...] = ef
        
        @property
        def ncoefsalpha(self):
            """
            Element ncoefsalpha ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 98
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__ncoefsalpha(self._handle)
        
        @ncoefsalpha.setter
        def ncoefsalpha(self, ncoefsalpha):
            _alegendre.f90wrap_bessel_expansion_data__set__ncoefsalpha(self._handle, \
                ncoefsalpha)
        
        @property
        def ncoefsalphap(self):
            """
            Element ncoefsalphap ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 98
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__ncoefsalphap(self._handle)
        
        @ncoefsalphap.setter
        def ncoefsalphap(self, ncoefsalphap):
            _alegendre.f90wrap_bessel_expansion_data__set__ncoefsalphap(self._handle, \
                ncoefsalphap)
        
        @property
        def coefsalpha(self):
            """
            Element coefsalpha ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 99
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__coefsalpha(self._handle)
            if array_handle in self._arrays:
                coefsalpha = self._arrays[array_handle]
            else:
                coefsalpha = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__coefsalpha)
                self._arrays[array_handle] = coefsalpha
            return coefsalpha
        
        @coefsalpha.setter
        def coefsalpha(self, coefsalpha):
            self.coefsalpha[...] = coefsalpha
        
        @property
        def coefsalphap(self):
            """
            Element coefsalphap ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 99
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__coefsalphap(self._handle)
            if array_handle in self._arrays:
                coefsalphap = self._arrays[array_handle]
            else:
                coefsalphap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__coefsalphap)
                self._arrays[array_handle] = coefsalphap
            return coefsalphap
        
        @coefsalphap.setter
        def coefsalphap(self, coefsalphap):
            self.coefsalphap[...] = coefsalphap
        
        @property
        def iptrsalpha(self):
            """
            Element iptrsalpha ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 101
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__iptrsalpha(self._handle)
            if array_handle in self._arrays:
                iptrsalpha = self._arrays[array_handle]
            else:
                iptrsalpha = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__iptrsalpha)
                self._arrays[array_handle] = iptrsalpha
            return iptrsalpha
        
        @iptrsalpha.setter
        def iptrsalpha(self, iptrsalpha):
            self.iptrsalpha[...] = iptrsalpha
        
        @property
        def iptrsalphap(self):
            """
            Element iptrsalphap ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 101
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__iptrsalphap(self._handle)
            if array_handle in self._arrays:
                iptrsalphap = self._arrays[array_handle]
            else:
                iptrsalphap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__iptrsalphap)
                self._arrays[array_handle] = iptrsalphap
            return iptrsalphap
        
        @iptrsalphap.setter
        def iptrsalphap(self, iptrsalphap):
            self.iptrsalphap[...] = iptrsalphap
        
        @property
        def ncoefsbeta1(self):
            """
            Element ncoefsbeta1 ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 102
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__ncoefsbeta1(self._handle)
        
        @ncoefsbeta1.setter
        def ncoefsbeta1(self, ncoefsbeta1):
            _alegendre.f90wrap_bessel_expansion_data__set__ncoefsbeta1(self._handle, \
                ncoefsbeta1)
        
        @property
        def ncoefsbeta2(self):
            """
            Element ncoefsbeta2 ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 102
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__ncoefsbeta2(self._handle)
        
        @ncoefsbeta2.setter
        def ncoefsbeta2(self, ncoefsbeta2):
            _alegendre.f90wrap_bessel_expansion_data__set__ncoefsbeta2(self._handle, \
                ncoefsbeta2)
        
        @property
        def coefsbeta1(self):
            """
            Element coefsbeta1 ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 103
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__coefsbeta1(self._handle)
            if array_handle in self._arrays:
                coefsbeta1 = self._arrays[array_handle]
            else:
                coefsbeta1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__coefsbeta1)
                self._arrays[array_handle] = coefsbeta1
            return coefsbeta1
        
        @coefsbeta1.setter
        def coefsbeta1(self, coefsbeta1):
            self.coefsbeta1[...] = coefsbeta1
        
        @property
        def coefsbeta2(self):
            """
            Element coefsbeta2 ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 103
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__coefsbeta2(self._handle)
            if array_handle in self._arrays:
                coefsbeta2 = self._arrays[array_handle]
            else:
                coefsbeta2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__coefsbeta2)
                self._arrays[array_handle] = coefsbeta2
            return coefsbeta2
        
        @coefsbeta2.setter
        def coefsbeta2(self, coefsbeta2):
            self.coefsbeta2[...] = coefsbeta2
        
        @property
        def iptrsbeta1(self):
            """
            Element iptrsbeta1 ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 104
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__iptrsbeta1(self._handle)
            if array_handle in self._arrays:
                iptrsbeta1 = self._arrays[array_handle]
            else:
                iptrsbeta1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__iptrsbeta1)
                self._arrays[array_handle] = iptrsbeta1
            return iptrsbeta1
        
        @iptrsbeta1.setter
        def iptrsbeta1(self, iptrsbeta1):
            self.iptrsbeta1[...] = iptrsbeta1
        
        @property
        def iptrsbeta2(self):
            """
            Element iptrsbeta2 ftype=integer pytype=int
            
            
            Defined at bessel_eval.fpp line 104
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _alegendre.f90wrap_bessel_expansion_data__array__iptrsbeta2(self._handle)
            if array_handle in self._arrays:
                iptrsbeta2 = self._arrays[array_handle]
            else:
                iptrsbeta2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _alegendre.f90wrap_bessel_expansion_data__array__iptrsbeta2)
                self._arrays[array_handle] = iptrsbeta2
            return iptrsbeta2
        
        @iptrsbeta2.setter
        def iptrsbeta2(self, iptrsbeta2):
            self.iptrsbeta2[...] = iptrsbeta2
        
        @property
        def dmemory(self):
            """
            Element dmemory ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 105
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__dmemory(self._handle)
        
        @dmemory.setter
        def dmemory(self, dmemory):
            _alegendre.f90wrap_bessel_expansion_data__set__dmemory(self._handle, dmemory)
        
        @property
        def time(self):
            """
            Element time ftype=double precision pytype=unknown
            
            
            Defined at bessel_eval.fpp line 105
            
            """
            return _alegendre.f90wrap_bessel_expansion_data__get__time(self._handle)
        
        @time.setter
        def time(self, time):
            _alegendre.f90wrap_bessel_expansion_data__set__time(self._handle, time)
        
        def __str__(self):
            ret = ['<bessel_expansion_data>{\n']
            ret.append('    ifbetas : ')
            ret.append(repr(self.ifbetas))
            ret.append(',\n    ifover : ')
            ret.append(repr(self.ifover))
            ret.append(',\n    ifsmall : ')
            ret.append(repr(self.ifsmall))
            ret.append(',\n    dnu1 : ')
            ret.append(repr(self.dnu1))
            ret.append(',\n    dnu2 : ')
            ret.append(repr(self.dnu2))
            ret.append(',\n    epsrequired : ')
            ret.append(repr(self.epsrequired))
            ret.append(',\n    k : ')
            ret.append(repr(self.k))
            ret.append(',\n    nintsab : ')
            ret.append(repr(self.nintsab))
            ret.append(',\n    nintscd : ')
            ret.append(repr(self.nintscd))
            ret.append(',\n    nintsef : ')
            ret.append(repr(self.nintsef))
            ret.append(',\n    ab : ')
            ret.append(repr(self.ab))
            ret.append(',\n    cd : ')
            ret.append(repr(self.cd))
            ret.append(',\n    ef : ')
            ret.append(repr(self.ef))
            ret.append(',\n    ncoefsalpha : ')
            ret.append(repr(self.ncoefsalpha))
            ret.append(',\n    ncoefsalphap : ')
            ret.append(repr(self.ncoefsalphap))
            ret.append(',\n    coefsalpha : ')
            ret.append(repr(self.coefsalpha))
            ret.append(',\n    coefsalphap : ')
            ret.append(repr(self.coefsalphap))
            ret.append(',\n    iptrsalpha : ')
            ret.append(repr(self.iptrsalpha))
            ret.append(',\n    iptrsalphap : ')
            ret.append(repr(self.iptrsalphap))
            ret.append(',\n    ncoefsbeta1 : ')
            ret.append(repr(self.ncoefsbeta1))
            ret.append(',\n    ncoefsbeta2 : ')
            ret.append(repr(self.ncoefsbeta2))
            ret.append(',\n    coefsbeta1 : ')
            ret.append(repr(self.coefsbeta1))
            ret.append(',\n    coefsbeta2 : ')
            ret.append(repr(self.coefsbeta2))
            ret.append(',\n    iptrsbeta1 : ')
            ret.append(repr(self.iptrsbeta1))
            ret.append(',\n    iptrsbeta2 : ')
            ret.append(repr(self.iptrsbeta2))
            ret.append(',\n    dmemory : ')
            ret.append(repr(self.dmemory))
            ret.append(',\n    time : ')
            ret.append(repr(self.time))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def bessel_eval(dnu, t):
        """
        alpha, alphader, vallogj, vallogy, valj, valy = bessel_eval(dnu, t)
        
        
        Defined at bessel_eval.fpp lines 127-271
        
        Parameters
        ----------
        dnu : unknown
        t : unknown
        
        Returns
        -------
        alpha : unknown
        alphader : unknown
        vallogj : unknown
        vallogy : unknown
        valj : unknown
        valy : unknown
        
        """
        alpha, alphader, vallogj, vallogy, valj, valy = \
            _alegendre.f90wrap_bessel_eval(dnu=dnu, t=t)
        return alpha, alphader, vallogj, vallogy, valj, valy
    
    @staticmethod
    def bessel_expeval(self, dnu, t):
        """
        aval, apval, bval1, bval2, valj, valy = bessel_expeval(self, dnu, t)
        
        
        Defined at bessel_eval.fpp lines 273-346
        
        Parameters
        ----------
        expdata : Bessel_Expansion_Data
        dnu : unknown
        t : unknown
        
        Returns
        -------
        aval : unknown
        apval : unknown
        bval1 : unknown
        bval2 : unknown
        valj : unknown
        valy : unknown
        
        """
        aval, apval, bval1, bval2, valj, valy = \
            _alegendre.f90wrap_bessel_expeval(expdata=self._handle, dnu=dnu, t=t)
        return aval, apval, bval1, bval2, valj, valy
    
    @staticmethod
    def bessel_evalabc(ifsmall, dnu, a, b, c):
        """
        bessel_evalabc(ifsmall, dnu, a, b, c)
        
        
        Defined at bessel_eval.fpp lines 348-362
        
        Parameters
        ----------
        ifsmall : int
        dnu : float
        a : float
        b : float
        c : float
        
        """
        _alegendre.f90wrap_bessel_evalabc(ifsmall=ifsmall, dnu=dnu, a=a, b=b, c=c)
    
    @staticmethod
    def bessel_phase_asymp(dnu, x, aval, apval):
        """
        bessel_phase_asymp(dnu, x, aval, apval)
        
        
        Defined at bessel_eval.fpp lines 364-395
        
        Parameters
        ----------
        dnu : float
        x : float
        aval : float
        apval : float
        
        """
        _alegendre.f90wrap_bessel_phase_asymp(dnu=dnu, x=x, aval=aval, apval=apval)
    
    @staticmethod
    def kummer_bessel_phase_asymp(dnu, t, aval, apval, appval):
        """
        kummer_bessel_phase_asymp(dnu, t, aval, apval, appval)
        
        
        Defined at bessel_eval.fpp lines 397-461
        
        Parameters
        ----------
        dnu : float
        t : float
        aval : float
        apval : float
        appval : float
        
        """
        _alegendre.f90wrap_kummer_bessel_phase_asymp(dnu=dnu, t=t, aval=aval, \
            apval=apval, appval=appval)
    
    @staticmethod
    def bessel_debye(dnu, t, bval1, bval2, valj, valy):
        """
        bessel_debye(dnu, t, bval1, bval2, valj, valy)
        
        
        Defined at bessel_eval.fpp lines 463-1022
        
        Parameters
        ----------
        dnu : float
        t : float
        bval1 : float
        bval2 : float
        valj : float
        valy : float
        
        """
        _alegendre.f90wrap_bessel_debye(dnu=dnu, t=t, bval1=bval1, bval2=bval2, \
            valj=valj, valy=valy)
    
    @staticmethod
    def bessel_taylor(dnu, t, alpha, alphader, vallogj, vallogy, valj, valy):
        """
        bessel_taylor(dnu, t, alpha, alphader, vallogj, vallogy, valj, valy)
        
        
        Defined at bessel_eval.fpp lines 1024-1077
        
        Parameters
        ----------
        dnu : float
        t : float
        alpha : float
        alphader : float
        vallogj : float
        vallogy : float
        valj : float
        valy : float
        
        """
        _alegendre.f90wrap_bessel_taylor(dnu=dnu, t=t, alpha=alpha, alphader=alphader, \
            vallogj=vallogj, vallogy=vallogy, valj=valj, valy=valy)
    
    @staticmethod
    def bessel_taylorj(dnu, t, val):
        """
        bessel_taylorj(dnu, t, val)
        
        
        Defined at bessel_eval.fpp lines 1079-1100
        
        Parameters
        ----------
        dnu : float
        t : float
        val : float
        
        """
        _alegendre.f90wrap_bessel_taylorj(dnu=dnu, t=t, val=val)
    
    @staticmethod
    def bessel_taylory(dnu, t, val):
        """
        bessel_taylory(dnu, t, val)
        
        
        Defined at bessel_eval.fpp lines 1102-1252
        
        Parameters
        ----------
        dnu : float
        t : float
        val : float
        
        """
        _alegendre.f90wrap_bessel_taylory(dnu=dnu, t=t, val=val)
    
    @staticmethod
    def bessel_taylorj_log(dnu, t, vallog):
        """
        bessel_taylorj_log(dnu, t, vallog)
        
        
        Defined at bessel_eval.fpp lines 1254-1282
        
        Parameters
        ----------
        dnu : float
        t : float
        vallog : float
        
        """
        _alegendre.f90wrap_bessel_taylorj_log(dnu=dnu, t=t, vallog=vallog)
    
    @staticmethod
    def bessel_taylorj_log2(dnu, t, dsign, vallog):
        """
        bessel_taylorj_log2(dnu, t, dsign, vallog)
        
        
        Defined at bessel_eval.fpp lines 1284-1317
        
        Parameters
        ----------
        dnu : float
        t : float
        dsign : float
        vallog : float
        
        """
        _alegendre.f90wrap_bessel_taylorj_log2(dnu=dnu, t=t, dsign=dsign, vallog=vallog)
    
    @staticmethod
    def bessel_taylory_log(dnu, t, vallog):
        """
        bessel_taylory_log(dnu, t, vallog)
        
        
        Defined at bessel_eval.fpp lines 1319-1452
        
        Parameters
        ----------
        dnu : float
        t : float
        vallog : float
        
        """
        _alegendre.f90wrap_bessel_taylory_log(dnu=dnu, t=t, vallog=vallog)
    
    @staticmethod
    def bessel_taylory_log0(dnu, t, vallog):
        """
        bessel_taylory_log0(dnu, t, vallog)
        
        
        Defined at bessel_eval.fpp lines 1454-1471
        
        Parameters
        ----------
        dnu : float
        t : float
        vallog : float
        
        """
        _alegendre.f90wrap_bessel_taylory_log0(dnu=dnu, t=t, vallog=vallog)
    
    @staticmethod
    def bessel_tensor_eval(ncoefs1, coefs1, ncoefs2, coefs2, iptr1, iptr2, a, b, c, \
        d, x, y, val1, val2):
        """
        bessel_tensor_eval(ncoefs1, coefs1, ncoefs2, coefs2, iptr1, iptr2, a, b, c, d, \
            x, y, val1, val2)
        
        
        Defined at bessel_eval.fpp lines 1473-1521
        
        Parameters
        ----------
        ncoefs1 : int
        coefs1 : unknown array
        ncoefs2 : int
        coefs2 : unknown array
        iptr1 : int
        iptr2 : int
        a : float
        b : float
        c : float
        d : float
        x : float
        y : float
        val1 : float
        val2 : float
        
        """
        _alegendre.f90wrap_bessel_tensor_eval(ncoefs1=ncoefs1, coefs1=coefs1, \
            ncoefs2=ncoefs2, coefs2=coefs2, iptr1=iptr1, iptr2=iptr2, a=a, b=b, c=c, \
            d=d, x=x, y=y, val1=val1, val2=val2)
    
    @staticmethod
    def bessel_chebs(x, n, pols):
        """
        bessel_chebs(x, n, pols)
        
        
        Defined at bessel_eval.fpp lines 1523-1581
        
        Parameters
        ----------
        x : unknown
        n : int
        pols : unknown array
        
        """
        _alegendre.f90wrap_bessel_chebs(x=x, n=n, pols=pols)
    
    @staticmethod
    def bessel_findint(nints, ab, x, int_bn, a, b):
        """
        bessel_findint(nints, ab, x, int_bn, a, b)
        
        
        Defined at bessel_eval.fpp lines 1583-1605
        
        Parameters
        ----------
        nints : int
        ab : unknown array
        x : unknown
        int_bn : int
        a : unknown
        b : unknown
        
        """
        _alegendre.f90wrap_bessel_findint(nints=nints, ab=ab, x=x, int_bn=int_bn, a=a, \
            b=b)
    
    @staticmethod
    def bessel_eval_init(dsize):
        """
        bessel_eval_init(dsize)
        
        
        Defined at bessel_eval.fpp lines 1607-1667
        
        Parameters
        ----------
        dsize : float
        
        """
        _alegendre.f90wrap_bessel_eval_init(dsize=dsize)
    
    @staticmethod
    def read_expansion(iw):
        """
        expdata = read_expansion(iw)
        
        
        Defined at bessel_eval.fpp lines 1669-1736
        
        Parameters
        ----------
        iw : int
        
        Returns
        -------
        expdata : Bessel_Expansion_Data
        
        """
        expdata = _alegendre.f90wrap_read_expansion(iw=iw)
        expdata = \
            f90wrap.runtime.lookup_class("alegendre.bessel_expansion_data").from_handle(expdata)
        return expdata
    
    @staticmethod
    def bessel_read_double_array_binary(iw, n, data):
        """
        bessel_read_double_array_binary(iw, n, data)
        
        
        Defined at bessel_eval.fpp lines 1738-1749
        
        Parameters
        ----------
        iw : int
        n : int
        data : unknown array
        
        """
        _alegendre.f90wrap_bessel_read_double_array_binary(iw=iw, n=n, data=data)
    
    @staticmethod
    def bessel_read_integer_array_binary(iw, n, idata):
        """
        bessel_read_integer_array_binary(iw, n, idata)
        
        
        Defined at bessel_eval.fpp lines 1751-1764
        
        Parameters
        ----------
        iw : int
        n : int
        idata : int array
        
        """
        _alegendre.f90wrap_bessel_read_integer_array_binary(iw=iw, n=n, idata=idata)
    
    @staticmethod
    def bessel_read_integer_binary(iw, idata):
        """
        bessel_read_integer_binary(iw, idata)
        
        
        Defined at bessel_eval.fpp lines 1766-1775
        
        Parameters
        ----------
        iw : int
        idata : int
        
        """
        _alegendre.f90wrap_bessel_read_integer_binary(iw=iw, idata=idata)
    
    @staticmethod
    def bessel_read_double_binary(iw, data):
        """
        bessel_read_double_binary(iw, data)
        
        
        Defined at bessel_eval.fpp lines 1777-1786
        
        Parameters
        ----------
        iw : int
        data : float
        
        """
        _alegendre.f90wrap_bessel_read_double_binary(iw=iw, data=data)
    
    @staticmethod
    def read_expansion16(iw):
        """
        expdata = read_expansion16(iw)
        
        
        Defined at bessel_eval.fpp lines 1788-1855
        
        Parameters
        ----------
        iw : int
        
        Returns
        -------
        expdata : Bessel_Expansion_Data
        
        """
        expdata = _alegendre.f90wrap_read_expansion16(iw=iw)
        expdata = \
            f90wrap.runtime.lookup_class("alegendre.bessel_expansion_data").from_handle(expdata)
        return expdata
    
    @staticmethod
    def bessel_read_double_array_binary16(iw, n, data):
        """
        bessel_read_double_array_binary16(iw, n, data)
        
        
        Defined at bessel_eval.fpp lines 1857-1868
        
        Parameters
        ----------
        iw : int
        n : int
        data : unknown array
        
        """
        _alegendre.f90wrap_bessel_read_double_array_binary16(iw=iw, n=n, data=data)
    
    @staticmethod
    def bessel_read_double_binary16(iw, data):
        """
        bessel_read_double_binary16(iw, data)
        
        
        Defined at bessel_eval.fpp lines 1870-1879
        
        Parameters
        ----------
        iw : int
        data : float
        
        """
        _alegendre.f90wrap_bessel_read_double_binary16(iw=iw, data=data)
    
    @staticmethod
    def bessel_gamma(x):
        """
        bessel_gamma = bessel_gamma(x)
        
        
        Defined at bessel_eval.fpp lines 1881-1904
        
        Parameters
        ----------
        x : float
        
        Returns
        -------
        bessel_gamma : unknown
        
        """
        bessel_gamma = _alegendre.f90wrap_bessel_gamma(x=x)
        return bessel_gamma
    
    _dt_array_initialisers = []
    

besseleval = Besseleval()

