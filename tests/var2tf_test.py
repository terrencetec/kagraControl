import sys
sys.path.append('../../')
from kagraControl.systemID import *
from control import *
def test_zpk(): # test tf2var with zpk entry and convert it back to tf with var2tf.
    z=([-1+1j,-1-1j,-3+5j,-3-5j,-0.01])
    p=([-3+1j,-3-1j,-5+6j,-5-6j,-10,-20,-30])
    k=3.14
    var_test,order_test,tf_test=tf2var(tf0=[z,p,k])
#     print(order_test)
#     print(tf_test.dcgain())
#     print(tf_test.pole())
#     print(tf_test.zero())
#     tf_test
    tf_test2=var2tf(var=var_test,order=order_test)
#     print(tf_test2.dcgain())
#     print(tf_test2.pole())
#     print(tf_test2.zero())
    print('The following tf printed should be 1/1')
    print(minreal(tf_test2/tf_test))

def test_numden(): # test tf2var with num/den entry and convert it back to tf with var2tf
    num=[1,2,3]
    den=[4,5,6,7]
    var_test,order_test,tf_test=tf2var(tf0=[num,den])
    # print(tf_test.dcgain())
    # print(tf_test.pole())
    # print(tf_test.zero())
    tf_test2=var2tf(var=var_test,order=order_test)
    # print(tf_test2.dcgain())
    # print(tf_test2.pole())
    # print(tf_test2.zero())
    print('The following tf printed should be 1/1')
    print(minreal(tf_test2/tf_test))

def test_var(): # test tf2var with var, order entry and convert it back to tf with var2tf
    var_test=[1,2,3,4,5,6,7,3.14]
    order_test=[1,2,2,2] # 1 simple zero, 1 simple pole, 2 complex zeros, 2 complex poles
    var_test2,order_test2,tf_test=tf2var(tf0=var_test,order=order_test)
    # print(tf_test.dcgain())
    # print(tf_test.pole())
    # print(tf_test.zero())
    tf_test2=var2tf(var=var_test2,order=order_test2)
    # print(tf_test2.dcgain())
    # print(tf_test2.pole())
    # print(tf_test2.zero())
    print('The following tf printed should be 1/1')
    print(minreal(tf_test2/tf_test))

print('Test tf2var(zpk)')
test_zpk()
print('Test tf2var(num/den)')
test_numden()
print('Test tf2var(var,order)')
test_var()
