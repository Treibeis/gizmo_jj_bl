import sys
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from txt import *
import random

select=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13])

x = np.array([250, 90, 31, 29.5, 10.49, 9.72, 8.33, 7.98, 3.41, 3.20, 1.42, 0.408, 0.045, 0.022])
y = np.array([2.725, 2.706, 2.573, 2.529, 2.738, 2.731, 2.743, 2.761, 2.770, 2.787, 3.271, 13.42, 3864, 20355])*1e3

C = np.array([[1]+[0]*13,[0, 370, 6.1,5.7,3.5,4.3,47.2,41.0,7.2,13]+[0]*4, [0,6.1,5800,2.2,1.7,1.9,13.4,11.8,2.7,4.1]+[0]*4,[0,5.7,2.2,2.4e4,1.6,1.8,12.6,11,2.5,3.9]+[0]*4,[0,3.5,1.7,1.6,30.3,1.19,6.16,5.48,1.76,2.5,1.6,53,8700,3.8e4],[0,4.3,1.9,1.8,1.19,30.1,7.72,6.79,2.04,2.86,1.7,55,9000,4e4],[0, 47.2,13.4,12.6,6.16,7.72,222,86,14.8,26.9,2.47,82,1.3e4,5.9e4],[0,41,11.8,11,5.48,6.79,86,173,13,23.4,2.1,71,1.1e4,5.1e4],[0,7.2,2.7,2.5,1.76,2.04,14.8,13,69,6.78,21,700,1.1e5,5e5],[0,13,4.1,3.9,2.5,2.86,26.9,23.4,6.78,93.9,24,780,1.3e5,5.6e5],[0]*4+[1.6,1.7,2.47,2.1,21,24,2.8e5,7300,1.2e6,5.3e6],[0]*4+[53,55,82,71,700,780,7300,1.2e7,3.9e7,8.4e9],[0]*4+[8700,9000,1.3e4,1.1e4,1.1e5,1.3e5,1.2e6,3.9e7,2.5e11,2.8e10],[0]*4+[3.8e4,4e4,5.9e4,5.1e4,5e5,5.6e5,5.3e6,8.4e9,2.8e10,2.7e13]])

yerr = np.array([C[i][i]**0.5 for i in range(len(y))])

para_default = [2.725,24.1e3,0.0,2.6]

def fit(mode, x=x[select], y=y[select], C=C[select][:,select], beta = -2.0, ref = 0.31, epochs = 50000, alpha = 1,para0=para_default,switch=True,nbatch=10,jp0=2000,fac=1):
	x_dim, y_dim = x.shape[0], y.shape[0]
	assert x_dim == y_dim
	y_tf = tf.placeholder(dtype=tf.float32)#, shape=[y_dim])
	x_tf = tf.placeholder(dtype=tf.float32)#, shape=[y_dim])
	C_tf = tf.placeholder(dtype=tf.float32)#, shape=[x_dim, y_dim])

	def clip_to_nonnegative(x):
		return tf.clip_by_value(x, 0., np.inf)
	def clip_to_ind(x):
		return tf.clip_by_value(x, 2.0, 3.0)
	ini0 = tf.constant_initializer(para0[0],dtype=tf.float32)
	ini1 = tf.constant_initializer(para0[1],dtype=tf.float32)
	ini2 = tf.constant_initializer(para0[2],dtype=tf.float32)
	inii = tf.constant_initializer(para0[3],dtype=tf.float32)
	dof = 14-int(switch)-(1-mode)-2
	sigma_chi = (2/dof)**0.5
	print('dof: {}, sigma_chi: {}'.format(dof, sigma_chi))
	with tf.variable_scope('test'+str(int(1e4*np.random.random()))) as scope:
		#scope.__init__(reuse=tf.AUTO_REUSE)
		T0 = tf.get_variable(name='T0'+str(switch),dtype=tf.float32,initializer=ini0,constraint=clip_to_nonnegative,shape=[])
		T1 = tf.get_variable(name='T1'+str(switch),dtype=tf.float32,initializer=ini1,constraint=clip_to_nonnegative,shape=[])
		if switch:
			T2 = tf.get_variable(name='T2'+str(switch),dtype=tf.float32,initializer=ini2,constraint=clip_to_nonnegative,shape=[])
		else:
			T2 = tf.get_variable(name='T2'+str(switch),dtype=tf.float32,initializer=ini2,constraint=clip_to_nonnegative,shape=[],trainable=False)
		if mode==0:
			index = tf.get_variable(name='index'+str(switch),dtype=tf.float32,initializer=inii,constraint=clip_to_ind,shape=[])
		else:
			index = tf.get_variable(name='index'+str(switch),dtype=tf.float32,initializer=inii,constraint=clip_to_ind,trainable=False,shape=[])

		T = T0 + T1*(x_tf/ref)**-index + T2*(x_tf/ref)**beta
		dis = T-y_tf
		dis_ = tf.tensordot(dis,tf.linalg.inv(C_tf),axes=[[0],[0]])
		obj = tf.tensordot(dis_,dis,axes=[[0],[0]])/(nbatch)#-2-int(mode==0)-int(switch==True))
	optimizer = tf.train.AdamOptimizer(alpha).minimize(obj)
	init_op = tf.global_variables_initializer()
	lT0, lT1, lT2, lindex = [], [], [], []
	lobj = []
	with tf.Session() as sess:
		sess.run(init_op)
		obj_val0 = 0
		for ep in range(epochs):
#			minibatch = [x for x in range(14)]
			minibatch = np.random.choice(x_dim, size=nbatch, replace=False)
			sess.run(optimizer, feed_dict={x_tf:x[minibatch], y_tf:y[minibatch], C_tf:C[minibatch][:,minibatch]})
			if ep%jp0==0:
				obj_val = sess.run(obj, feed_dict={x_tf:x, y_tf:y, C_tf:C})
				dis_val = sess.run(dis, feed_dict={x_tf:x, y_tf:y, C_tf:C})
				dis_val_ = np.dot(np.linalg.inv(C),dis_val)#sess.run(dis_, feed_dict={x_tf:x, y_tf:y, C_tf:C})
				print('epoch {} obj = {}'.format(ep, 14*obj_val))
				if abs(obj_val-obj_val0)/obj_val<1e-6:
					break
				else:
					obj_val0 = obj_val
					lT0.append(sess.run(T0))
					lT1.append(sess.run(T1))
					lT2.append(sess.run(T2))
					lindex.append(sess.run(index))
					lobj.append(obj_val)
		obj_val = sess.run(obj, feed_dict={x_tf:x, y_tf:y, C_tf:C})
		dis_val = sess.run(dis, feed_dict={x_tf:x, y_tf:y, C_tf:C})
		T_pre = sess.run(T, feed_dict={x_tf:x, y_tf:y, C_tf:C})
		print('epoch {} obj = {}'.format(ep, 14*obj_val))
		#print(dis_val)
	lobj = np.array(lobj)
	o = lobj< (obj_val*14/dof + sigma_chi)*dof/14
	lT0 = np.array(lT0)
	lT1 = np.array(lT1)
	lT2 = np.array(lT2)
	lindex = np.array(lindex)
	return [np.mean(lT0[o]), np.mean(lT1[o]), np.mean(lT2[o]), np.mean(lindex[o]), np.ptp(lT0[o])/2, np.ptp(lT1[o])/2, np.ptp(lT2[o])/2, np.ptp(lindex[o])/2]# np.std(lT0[o]), np.std(lT1[o]), np.std(lT2[o]), np.std(lindex[o])]

if __name__ == '__main__':
	#print(x.shape,y.shape,C.shape)
	#print(C[:5,:5])
	tag = 0
	ref = 0.31
	beta = -2.0
	mode = int(sys.argv[1])
	nbatch = int(sys.argv[2])
	#print(np.linalg.inv(C))
	#print(np.dot(np.linalg.inv(C),C))
	
	if tag==0:
		out0 = fit(mode,alpha=1,epochs=100000,para0=[2.725,24.1e3,0.0,2.599],ref=ref,beta=beta,switch=True,nbatch=nbatch)
		out1 = fit(mode,alpha=1,epochs=100000,para0=[2.725,24.1e3,0.0,2.599],ref=ref,beta=beta,switch=False,nbatch=nbatch)
		print('Without free-free, [T0, T1, T2, index]: {} [mK]'.format(out1))
		print('With free-free, [T0, T1, T2, index]: {} [mK]'.format(out0))
		totxt('fit_'+str(mode)+'.txt',[out0+[ref,beta],out1+[ref,beta]],0,0,0)
	else:
		out = retxt('fit_'+str(mode)+'.txt',2,0,0)
		out0 = out[0]
		out1 = out[1]
		ref, beta = out[0][8], out[0][9]
	
	lx = np.logspace(-2,np.log10(300),1000)
	def model(x, para, ref = ref, beta = beta):
		T0, T1, T2, index = para[0], para[1], para[2], para[3]
		return T0 + T1*(x/ref)**-index + T2*(x/ref)**beta

	plt.figure()
	plt.plot(lx[lx>2],model(lx, out1)[lx>2],label='Without free-free')
	plt.plot(lx[lx>2],model(lx, out0)[lx>2],'--',label='With free-free: '+r'$T_{\mathrm{ff}}('+str(ref)+'\ [\mathrm{GHz}])='+str(int(out0[2]))+'\ \mathrm{mk}$')
	plt.errorbar(x[x>3],y[x>3],yerr=yerr[x>3],label='Data',marker='o')
	plt.xlabel(r'$\nu\ [\mathrm{GHz}]$')
	plt.ylabel(r'$T\ [\mathrm{mK}]$')
	plt.xlim(2,300)
	plt.xscale('log')
	#plt.yscale('log')
	plt.legend()
	plt.tight_layout()
	plt.savefig('fit0_'+str(mode)+'.pdf')

	plt.figure()
	plt.plot(lx[lx<3],model(lx, out1)[lx<3],label='Without free-free')
	plt.plot(lx[lx<3],model(lx, out0)[lx<3],'--',label='With free-free: '+r'$T_{\mathrm{ff}}('+str(ref)+'\ [\mathrm{GHz}])='+str(int(out0[2]))+'\ \mathrm{mk}$')
	plt.errorbar(x[x<3],y[x<3],yerr=yerr[x<3],label='Data',marker='o')
	plt.xlim(1e-2,3)
	plt.xlabel(r'$\nu\ [\mathrm{GHz}]$')
	plt.ylabel(r'$T\ [\mathrm{mK}]$')
	plt.xscale('log')
	plt.yscale('log')
	plt.legend(prop={'size':12})
	plt.tight_layout()
	plt.savefig('fit1_'+str(mode)+'.pdf')

	plt.figure()
	plt.plot(x,(model(x, out1)-y)/yerr,label='Without free-free')
	plt.plot(x,(model(x, out0)-y)/yerr,'--',label='With free-free: '+r'$T_{\mathrm{ff}}('+str(ref)+'\ [\mathrm{GHz}])='+str(int(out0[2]))+'\ \mathrm{mk}$')
	plt.plot(lx,[0 for i in lx],'k-.',lw=1)
	plt.xlim(1e-2,300)
	plt.xlabel(r'$\nu\ [\mathrm{GHz}]$')
	plt.ylabel(r'$\Delta T/\sigma_{T}$')
	plt.legend()
	plt.xscale('log')
	plt.tight_layout()
	plt.savefig('error_'+str(mode)+'.pdf')

	print(out0[2]*(3.15/ref)**-2 / (out0[1]*(3.15/ref)**-out0[3]))
	plt.show()
	
	
