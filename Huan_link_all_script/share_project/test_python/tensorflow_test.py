import tensorflow as tf
# 此处代码需要使用 tf 1 版本运行
a_ph = tf.placeholder(tf.float32, name='variable_a')
b_ph = tf.placeholder(tf.float32, name='variable_b')

c_op = tf.add(a_ph, b_ph, name='variable_c')

sess = tf.InteractiveSession()
init = tf.global_variables_initializer()
sess.run(init)
c_numpy = sess.run(c_op, feed_dict={a_ph: 2., b_ph: 4.})
print('a+b=',c_numpy)

# 此处代码需要使用 tf 2 版本运行
a = tf.constant(2.)
b = tf.constant(4.)

print('a+b=', a+b)