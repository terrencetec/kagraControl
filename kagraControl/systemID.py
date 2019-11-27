import json
import numpy as np
from control import *
from scipy.optimize import *
# sample_complex=[np.pi+np.pi*1j]*10
# sample_f=np.linspace(0.01,10,10).tolist()

class measurement:
    # TODO make it work with .xml dtt measurements
    ''' A measurement objects which store the type of measurements.
        Usage:  Specify a json file which has all the required keys and items. measurement(filename='myMeasurement.json')
                Specify a dictionary which has all the required keys and items. measurement(dict=myDict), myDict={'type'='ts',input=[...],output=[...]}
                Specify the type of measurement: time series 'ts' or 'time series', complex numbers 'cplx' or 'complex' or frequency response 'fr' or 'frequency response'. measurement(type='ts',...)

                If type is time series, then the time axis t and input/output data (of the process) should be specified as well.
                if input data is not specified, a while noise of amplitude 1 will be assumed (?). measurement(type='ts',t=myTimeAxis,input=myInputData,output=myOutputData)
                If measurement data is composed of multiple averages, then input/output should be a list of matching data.

                If type is complex, then the frequency axis and complex numbers should be specfied. measurement(type='cplx',f=myFreqAxis,tf_complex=myComplexNumbers)
                Supplying the real parts and the imaginary parts of the complex numbers separately is also supported. measurement(type='cplx',f=myFreqAxis,tf_real=myRealNumbers,tf_imag=myImagNumbers)

                If type is fr, then frequency axis and mag/phase should be specfied. measurement(type='fr',f=myFreqAxis,mag=myMagData,phase=myPhaseData)
    '''
    def __init__(self,**kwargs):
        # TODO: check if the data has the same len as the axis.
        try:                            #Try of the keyword 'type' is used. if yes, then proceed to unpack
            self.type=kwargs['type']
            if self.type == 'ts' or self.type == 'time series':       #
                try:
                    self.t=kwargs['t']
                except:
                    print('Missing Expected argument: t.')
                    return(None)
                try:
                    self.output=kwargs['output']
                except:
                    print('Missing Expected argument: output.')
                    return(None)
                if len(np.shape(self.output))==1:
                    self.output=np.array([self.output])
                self.averages=np.shape(self.output)[0] # Maybe not useful.
                try:
                    self.input=kwargs['input']
                except:
                    print('Missing Expected argument: input. Proceed with default value 1.')
                    self.input=np.ones_like(self.output) # Check if .tolist() is very useful. Checked, it is required for json to work. # .tolist() is done when doing json save.
                if len(np.shape(self.input))==1:
                    self.input=np.array([self.input])
            elif self.type == 'cplx' or self.type == 'complex':
                try:
                    # if type(kwargs['f'])==np.ndarray:
                    #     self.f=kwargs['f'].tolist()
                    # else:
                    self.f=kwargs['f']
                except:
                    print('Missing Expected argument: f.')
                    raise
                # try:
                #     # self.tf_real=np.real(kwargs['tf_complex']).tolist()
                #     # self.tf_imag=np.imag(kwargs['tf_complex']).tolist()
                # except:
                try:
                    self.tf_real=kwargs['tf_real']
                    self.tf_imag=kwargs['tf_imag']
                    if len(np.shape(self.tf_real))==1:
                        self.tf_real=np.array([self.tf_real])
                    if len(np.shape(self.tf_imag))==1:
                        self.tf_imag=np.array([self.tf_imag])
                    self.averages=np.shape(self.tf_real)[0] # Maybe not useful.
                    if self.f[0]==0:
                        print('Zero Frequency data omitted')
                        self.f=self.f[1:len(self.f)]
                        if type(self.tf_real)==np.ndarray:
                            self.tf_real=self.tf_real.tolist()
                        if type(self.tf_imag)==np.ndarray:
                            self.tf_imag=self.tf_imag.tolist()
                        for i in range(len(self.tf_real)):
                            self.tf_real[i]=self.tf_real[i][1:len(self.tf_real[i])]
                        for i in range(len(self.tf_imag)):
                            self.tf_imag[i]=self.tf_imag[i][1:len(self.tf_imag[i])]
                except:
                    print('Missing Expected argument: tf_complex or tf_real and tf_imag.')
                    raise
            elif self.type == 'fr' or self.type == 'frequency response':
                print('Type fr not supported yet, try using time series or complex numbers')
                return(None)
        except: # if argument is a json filename or dictionary.
            try:
                self.__dict__=read_measurement(filename=kwargs['filename'])
                # unpack_dict(self)
            except:
                try:
                    self.__dict__=kwargs['dict']
                    # unpack_dict(kwargs)
                except:
                    print('No data in the correct format given. Proceed with an empty measurement object.')
        try:
            self.bounds
        except:
            self.bounds=None
        try:
            self.tf0
        except:
            self.tf0=None
        try:
            self.order
        except:
            self.order=None
        try:
            self.min_method
        except:
            self.min_method=None
        try:
            self.prefilter
        except:
            self.prefilter=[[1],[1]]

    def zpkfit(
            self,
            minimizer_args={}
            ):
        self.var_fit,self.order,_=fit_tf(measObj=self,tf0=self.tf0,order=self.order,bounds=self.bounds,min_method=self.min_method,prefilter=self.prefilter,minimizer_args=minimizer_args) # implement callback for minimzer
        if type(self.var_fit)==np.ndarray:
            self.var_fit=self.var_fit.tolist()


def tf2var(tf0,order=None):
    if type(tf0)==TransferFunction:
        # Extract critical frequencies, Q-factors and gain.
        mask_simple_z=tf0.zero().imag==0
        mask_complex_z=tf0.zero().imag!=0
        mask_simple_p=tf0.pole().imag==0
        mask_complex_p=tf0.pole().imag!=0
        simple_z_f=(np.abs(tf0.zero()[mask_simple_z])/2/np.pi).tolist()
        simple_p_f=(np.abs(tf0.pole()[mask_simple_p])/2/np.pi).tolist()
        complex_z_f=(np.abs(tf0.zero()[mask_complex_z])/2/np.pi).tolist()
        complex_z_Q=(-1/2/((tf0.zero()[mask_complex_z]).real/np.abs(tf0.zero()[mask_complex_z]))).tolist()
        complex_p_f=(np.abs(tf0.pole()[mask_complex_p])/2/np.pi).tolist()
        complex_p_Q=(-1/2/((tf0.pole()[mask_complex_p]).real/np.abs(tf0.pole()[mask_complex_p]))).tolist()
        gain=tf0.dcgain().tolist()
        var=simple_z_f+simple_p_f
        for i in range(len(complex_z_f)):
            if (i%2)==0: # even number to prevent double counting.
                var+=[complex_z_f[i]]
                var+=[complex_z_Q[i]]
        for i in range(len(complex_p_f)):
            if (i%2)==0: # even number
                var+=[complex_p_f[i]]
                var+=[complex_p_Q[i]]
        var+=[gain]
        order=[len(simple_z_f),len(simple_p_f),int(len(complex_z_f)),int(len(complex_p_f))]
        return(var,order,tf0)
    elif type(tf0)==tuple or type(tf0)==list:
        # print(tf0)

        if np.shape(tf0)[0]==2: #num/den
            num=tf0[0]
            den=tf0[1]
            tf1=tf(num,den)
            return(tf2var(tf1))
        elif np.shape(tf0)[0]==3: # zpk FIXME: deal with imaginary terms so they go in pairs.
            z=tf0[0]
            p=tf0[1]
            k=tf0[2]
            tf1=tf([1],[1])
            for i in range(len(z)):
                tf1*=tf([1,-z[i]],[1])
            for i in range(len(p)):
                tf1*=tf([1],[1,-p[i]])
            tf1*=tf([k],[1])/tf1.dcgain().tolist()
            return(tf2var(tf1))
        elif order != None: # np.shape(tf0)[0]==1: # tf0 is in the form of [f1,f2...,fn,Qn,k] and order is given.
            var1=tf0
            tf1=var2tf(var=var1,order=order)
            # print(tf1)
            return(tf2var(tf1))

def var2tf(var,order):
    tf0=tf([1],[1])
    skipOne=False
    for i in range(len(var)-1):
        if i<order[0]:
            tf0*=tf([1,var[i]*2*np.pi],[1])
        elif i<sum(order[0:2]):
            tf0*=tf([1],[1,var[i]*2*np.pi])
        elif i<sum(order[0:3]):
            if skipOne==True:
                skipOne=False
            else:
                fn=var[i]
                Q=var[i+1]
                tf0*=tf([1,(2*np.pi*fn)/Q,(2*np.pi*fn)**2],[(2*np.pi*fn)**2])
                skipOne=True
        else:
            if skipOne==True:
                skipOne=False
            else:
                fn=var[i]
                Q=var[i+1]
                tf0*=tf([(2*np.pi*fn)**2],[1,(2*np.pi*fn)/Q,(2*np.pi*fn)**2])
                skipOne=True
    tf0*=var[-1]/tf0.dcgain().tolist()
    # print('hi')
    return(tf0)
default={
        'min_Q':0.5,
        'max_Q':1e6,
        'min_k':-1e3,
        'max_k':1e3
        }
# fit measurement with a transfer function, given initial guess and bounds. If initial guess is not given, global methods will be used which is slow. In this case, either order or bounds must be specified.
def fit_tf(measObj, # Non-empty systemID.measurement object
        tf0=None, # Initial guess. control.tf object (SISO only), or 2-tuple (num/den) or 3-tuple(zpk), or list [f1,f2,...fn,Qn,K] work with the order argument
        order=None, # list of shape(4,), [No. of simple zeros, No. of simple poles, No. of complex zeros,  No. of complex poles]. Useful when tf0 is not specfied.
        bounds=None, # list of shape(2,) lower and upper bound of all critical frequencies [Hz] and Q values. list of shape (n,2), n has to match the total number of poles, zeros..
        min_method=None, # methods used in scipy.optimize.minimize or alternatively can be used to call global optimization from scipy.
        prefilter=[[1],[1]], # filter apply to the measured data before fitting, e.g. whitening. specification same as tf0.
        minimizer_args={}):
    prefilter_var,prefilter_order,_=tf2var(prefilter)
    prefilter=var2tf(prefilter_var,prefilter_order)
    # print(prefilter)
    if tf0!=None: # convert tf0 to variable form.
        tf0_var,_,__=tf2var(tf0,order)
        # tf0=var2tf(tf0_var,order)
    # print(tf0_var)
        minimizer_args['x0']=tf0_var
    # if tf0==None: # No initial guess specfied, use order, f_bounds and global methods.
    if order==None: # Terminate because those must be specified.
        print('If initial guess tf0 is not specified, order must be specified.')
        return(None)
    else:
        if bounds==None:
            print('bounds not specified, using measurement frequency bounds.')
            bounds=[[min(measObj.f),max(measObj.f)]]*order[0]
            bounds+=[[min(measObj.f),max(measObj.f)]]*order[1]
            skipOne=False
            for _ in range(order[2]+order[3]):
                if skipOne==True:
                    skipOne=False
                else:
                    bounds+=[[min(measObj.f),max(measObj.f)]]
                    bounds+=[[default['min_Q'],default['max_Q']]]
                    skipOne=True
            bounds+=[[default['min_k'],default['max_k']]]
        ## unpack global optimizer arguments here.
        if measObj.type in ['ts','time series']:
            t=measObj.t
            input=measObj.input
            output=measObj.output
            res=dual_annealing(cost_fit_ts,bounds=bounds,args=(t,input,output,order,prefilter),**minimizer_args)
            return(res.x,order,var2tf(var=res.x,order=order))
        elif measObj.type in ['cplx','complex']:
            f=measObj.f
            tf_real=measObj.tf_real
            tf_imag=measObj.tf_imag
            tf_real=np.array(tf_real)
            tf_imag=np.array(tf_imag)
            tf_complex=tf_real+1j*tf_imag
            # print(bounds)
            res=dual_annealing(cost_fit_cplx,bounds=bounds,args=(f,tf_complex,order,prefilter),**minimizer_args)
            return(res.x,order,var2tf(var=res.x,order=order))
        elif measObj.type in ['fr', 'frequency response']: # Convert to complex numbers and do it like complex.
            print('Type fr not supported yet, try using time series or complex numbers')
            return(None)
        else:
            print('Please specify a systemID.measurement object for fitting.')
            return(None)

# cost function for fitting
def cost_fit_ts(
                var, # variables to be adjusted
                t,
                input,
                output,
                order,
                prefilter=1
                ):
    t=np.linspace(t[0],t[-1],len(t))
    tf0=var2tf(var=var,order=order)*prefilter # FIXME: make sure prefilter is a transfer function object.
    if type(output)!=np.ndarray:
        output=np.array(output)
    residual=0
    for i in range(len(output)):
        _,output_fit,__=forced_response(sys=tf0,T=t,U=input[i])
        residual=np.sum(((output[i]-output_fit)/output[i])**2)
    return(residual)

def cost_fit_cplx(
                var, # variables to be adjusted
                f,
                tf_complex,
                order,
                prefilter=1
                ):
    # t=np.linspace(t[0],t[-1],len(t))
    # print(var,order)
    tf0=var2tf(var=var,order=order)*prefilter # FIXME: make sure prefilter is a transfer function object.
    if type(f)!=np.ndarray:
        f=np.array(f)
    if type(tf_complex)!=np.ndarray:
        tf_complex=np.array(tf_complex)
    residual=0
    for i in range(len(tf_complex)):
        tf_complex_fit=tf0.horner(1j*np.pi*2*f)[0][0]
        # print(np.abs((tf_complex[i]-tf_complex_fit)/np.abs(tf_complex[i])))
        residual=np.sum(np.abs((tf_complex[i]-tf_complex_fit)/np.abs(tf_complex[i])))
    return(residual)
    # TODO: consider putting these below to utilities and import them instead. So the whole problem is not excessively long.
def read_measurement(filename):
    '''Read a json file and return which can be unpacked to self.type, self.t,... etc'''
    if '.json' not in filename:
        filename=filename+'.json'
    with open(filename) as file:
        varDict=json.load(file)
    return(varDict)

def write_measurement(measObj,filename):
    '''Write self.__dict__ to a json file according to the measurement object specfied. Returns nothing'''
    if '.json' not in filename:
        filename=filename+'.json'
    for key in measObj.__dict__.keys():
        if type(measObj.__dict__[key])==np.ndarray:
            measObj.__dict__[key]=measObj.__dict__[key].tolist()
    # pack_dict(measObj)
    with open(filename,'w') as file:
        json.dump(measObj.__dict__,file)

# def unpack_dict(measObj):
#     '''unpack the dictionary to self.type, self.t... in the measurement object.'''
#     kwargs=measObj.dict
#     measObj.type=kwargs['type']
#     if measObj.type == 'ts' or measObj.type == 'time series':       #
#         try:
#             measObj.t=kwargs['t']
#         except:
#             print('Missing Expected argument: t.')
#             return(None)
#         try:
#             measObj.output=kwargs['output']
#         except:
#             print('Missing Expected argument: output.')
#             return(None)
#         if len(np.shape(measObj.output))==1:
#             measObj.output=[measObj.output]
#         measObj.averages=np.shape(measObj.output)[0] # Maybe not useful.
#         try:
#             measObj.input=kwargs['input']
#         except:
#             print('Missing Expected argument: input. Proceed with default value 1.')
#             measObj.input=np.ones_like(measObj.output).tolist() # Check if .tolist() is very useful.
#     elif measObj.type == 'cplx' or measObj.type == 'complex':
#         try:
#             measObj.f=kwargs['f']
#         except:
#             print('Missing Expected argument: f.')
#             return(None)
#         try:
#             measObj.tf_complex=kwargs['tf_complex']
#         except:
#             print('Missing Expected argument: tf_complex.')
#             return(None)
#         if len(np.shape(measObj.tf_complex))==1:
#             measObj.tf_complex=[measObj.tf_complex]
#         measObj.averages=np.shape(measObj.tf_complex)[0] # Maybe not useful.
#     elif measObj.type == 'fr' or measObj.type == 'frequency response':
#         print('Type fr not supported yet, try using time series or complex numbers')
#         return(None)
#
# def pack_dict(measObj):
#     measObj.dict=measObj.__dict__
    # measObj.dict['type']=measObj.type
    # if measObj.type in ['ts','time series']:
    #     measObj.dict['t']=measObj.t
    #     measObj.dict['input']=measObj.input
    #     measObj.dict['output']=measObj.output
    # elif measObj.type in ['cplx','complex']:
    #     measObj.dict['f']=measObj.f
    #     measObj.dict['tf_complex']=measObj.tf_complex
    # elif measObj.type in ['fr','frequency response']:
    #     measObj.dict['f']=measObj.f
    #     measObj.dict['mag']=measObj.mag
    #     measObj.dict['phase']=measObj.phase
    # measObj.dict['']
def test_systemID():
    return('systemID function test')

# sample_measurement=measurement(type='cplx',tf_complex=sample_complex,f=sample_f)
