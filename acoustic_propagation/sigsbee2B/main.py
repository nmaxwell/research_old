
from setup import *
temp_dir = "/workspace/tmp/research/"

(n1, n2) = shape(model)

session = {}
session['grid'] = grid


G = lambda x : exp(-1.0/x) if x > 0.0 else 0.0

def damping(x,y, grid):
    y1 = 0.0
    y2 = 8950.0
    x1 = 3700.0
    x2 = 27000.0
    
    ys = 50.0
    xs = 50.0
    
    return 8000.0*( G((y1-y)/ys) + G((y-y2)/ys) + G((x1-x)/xs) + G((x-x2)/xs) )

"""
import pylab as p
Y = grid.X2[0]
x = grid.X1[400][0]
Z = model[400]

damp = [ damping(x,y,grid) for y in Y ]

p.plot(Y,damp)
p.plot(Y,Z)
p.show()
quit()
"""

"""
import pylab as p
X = transpose(grid.X1)[200]
y = (grid(0,200))[1]
Z = transpose(model)[200]

damp = [ damping(x,y,grid) for x in X ]

p.plot(X,damp)
p.plot(X,Z)
p.show()
quit()
"""







session['drive_index1'] = session['grid'].index1( 1600 )
session['drive_index2'] = session['grid'].index2( 1000 )

session['initial_condition'] = lambda x,y: exp(-norm( array((x,y)) - array((16000.0, 8300.0)) )**2/(2.0*100.0**2))

session['driving_function'] = lambda t: 0.0
session['velocity_function'] = lambda x,y: 0.0
session['damping_function'] = damping
session['velocity'] = model
session['damping'] = array([])
session['driving'] = []
session['expansion_order'] = 3
session['hdaf_order'] = 8
session['hdaf_gamma'] = 0.8
session['time_step'] = 0.01
session['final_time'] = 1.0
session['u_fnames'] = []
session['time'] = []

P = propagator1()

def init(session, P ):
    print 'initializing, evaluating damping'
    session['damping'] = session['grid'].evaluate( lambda x,y: (session['damping_function'])(x,y,session['grid']) )
    print 'initializing propagator'
    P.set( session['grid'], velocity=session['velocity'], damping=session['damping'],  expansion_order=session['expansion_order'],   hdaf_order=session['hdaf_order'], hdaf_gamma=session['hdaf_gamma'] )
    write_png(session['damping'], "damp.png", major_scale=std(session['damping'])*3, center=mean(session['damping']), ordering='rm' )

def run( session ):
    
    grid = session['grid']
    
    u = grid.evaluate( session['initial_condition'] )
    v = grid.zeros()
    driving_1 = grid.zeros()
    driving_2 = grid.zeros()
    i0 = session['drive_index1']
    j0 = session['drive_index2']
    print i0, j0
    fname = temp_dir  + str(time_module.time())
    scipy.io.savemat( fname, {'m':u} )
    #results[str(step)] = { 'time': time, 'step':step, 'fname':fname }
    u_fnames = [fname]
    driving = []
    time_step = session['time_step']
    final_time = session['final_time']
    step=0
    time = [step]
    
    while time[step] <= final_time:
        
        print "step:", step, "\ttime:", time[step], "\twall time", wall_time(), "\tmagnitude:", norm(u)
        
        driving_1[i0,j0] = (session['driving_function'])(time[step])
        driving_2[i0,j0] = (session['driving_function'])(time[step]+time_step)
        driving.append(driving_1[i0,j0])
        
        P( time_step, u,v, u,v, driving_1, driving_2 )
        
        fname = temp_dir  + str(time_module.time())
        scipy.io.savemat( fname, {'m':u} )
        u_fnames.append(fname)
        
        step += 1
        time.append(time_step*step)
        
        session['u_fnames'] = u_fnames
        session['driving'] = driving
        session['time'] = time

    


def FT( session, frequency ):
    
    grid = session['grid']
    time_step = session['time_step']
    FT_sum_real = grid.zeros()  
    FT_sum_imag = grid.zeros()
    
    for step,time in enumerate(session['time']):
        print time
        data = scipy.io.loadmat( session['u_fnames'][step] )
        u = data['m']
        
        FT_sum_real += u*cos(-2.0*pi*time*frequency)*time_step
        FT_sum_imag += u*sin(-2.0*pi*time*frequency)*time_step
    
    return (FT_sum_real, FT_sum_imag)



def save( session, fname ):
    del(session['velocity_function'])
    del(session['damping_function'])
    del(session['driving_function'])
    del(session['initial_condition'])
    output = open(fname, 'wb')
    pickle.dump(session, output)
    output.close()


def load( fname ):
    input = open(fname, 'rb')
    session = pickle.load(input)
    input.close()
    return session

def render( data, fname ):
    ar = data.copy()
    write_png( ar, fname, center=mean(ar), major_scale=std(ar)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )

def render_all( session, directory="/workspace/output/acoustic_propagation/" ):
    
    for step,time in enumerate(session['time']):
        print time
        data = scipy.io.loadmat( session['u_fnames'][step] )
        u = data['m']
        
        render( u, directory+"%04d.png"%step )

if __name__ == "__main__":
    
    print "running main"
    
    init(session, P)
    run(session)
    save(session, 'session')
    
    render_all(session)

