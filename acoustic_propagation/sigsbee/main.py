
from setup import *

dir = "/tmp/research/"



def ricker_wavelet(t, sigma, gamma, tau ):
    return -sqrt(2.0/pi)*sigma*gamma*(sigma-2.0*sigma*gamma*(sigma*t-tau)**2)*exp(-gamma*(sigma*t-tau)**2)



(n1, n2) = shape(model)

session = {}
session['grid'] = grid2d( dx1=1.0, dx2=1.0, n1=n1, n2=n2)
grid = session['grid']



def damping(x,y, grid):
    scale = 50.0
    return 10000*( cosh(abs(x-grid.a1)/scale)**-2 + cosh(abs(x-grid.b1)/scale)**-2 + cosh(abs(y-grid.a2)/scale)**-2 + cosh(abs(y-grid.b2)/scale)**-2  )

session['drive_index1'] = session['grid'].index1( 1000 )
session['drive_index2'] = session['grid'].index2( 1000 )
session['driving_function'] = lambda t: exp(-(t-0.1)**2*5.0)
#ricker_wavelet(t=t, sigma=1.5*20, gamma=8, tau=1 )
session['velocity_function'] = lambda x,y: velocity((x,y))
session['damping_function'] = damping
session['velocity'] = model
session['damping'] = array([])
session['driving'] = []
session['expansion_order'] = 2
session['m1'] = 8
session['m2'] = 8
session['gamma1'] = 0.8
session['gamma2'] = 0.8
session['time_step'] = 0.005
session['final_time'] = 1.0
session['u_fnames'] = []
session['time'] = []

P = propagator1()

def init(session, P ):
    print 'initializing, evaluating damping'
    session['damping'] = session['grid'].evaluate( lambda x,y: (session['damping_function'])(x,y,session['grid']) )
    print 'initializing propagator'
    P.set( session['grid'], velocity=session['velocity'], damping=session['damping'],  expansion_order=session['expansion_order'],   hdaf_m1=session['m1'], hdaf_m2=session['m2'], hdaf_gamma1=session['gamma1'], hdaf_gamma2=session['gamma2'] )


def run( session ):
    
    grid = session['grid']
    
    u = grid.zeros()
    v = grid.zeros()
    driving_1 = grid.zeros()
    driving_2 = grid.zeros()
    i0 = session['drive_index1']
    j0 = session['drive_index2']
    print i0, j0
    fname = dir  + str(time_module.time())
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
        
        fname = dir  + str(time_module.time())
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

def render_all( session, dir="/workspace/output/scratch/" ):
    
    for step,time in enumerate(session['time']):
        print time
        data = scipy.io.loadmat( session['u_fnames'][step] )
        u = data['m']
        
        render( u, dir+"%04d.png"%step )

if __name__ == "__main__":
    
    print "running main"
    
    init(session, P)
    init(session, P)
    run(session)
    save(session, 'session')
    
    
    
    render_all(session)



