import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from scipy import stats

#planets = ['Sun','Earth']
planets = ['Sun','Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter','Saturn','Uranus','Neptune','Pluto']
moons = []
#moons = ['Ganymede','Titan','Callisto','Io','Luna','Europa','Triton','Titania','Rhea','Oberon']
body_list = planets + moons
number_bodies = len(body_list)

#==============================================================================
# plot solar system
#==============================================================================

#for body in body_list:
#    x, y = np.loadtxt('verlet_'+body+'.txt', delimiter=',', unpack=True)
#    plt.plot(x,y, label=body)

#dt = '15 minutes'
#time_span = '250 years'

#plt.xlabel('x [AU]')
#plt.ylabel('y [AU]')
#plt.title('Solar system, time step = {}, time span = {}'.format(dt,time_span))
#plt.legend(body_list)
#plt.show(block=False)

#==============================================================================
# plot x,y earth-sun system
#==============================================================================

#dts = ['1 month','1 week', '1 day', '1 hour', '15 minutes', '1 minute', '1 s']
#for i in range(0,5):
#    plt.figure()
#    x, y = np.loadtxt('verlet_dt_equals_'+str(i)+'.txt', delimiter=',', unpack=True)
#    plt.plot(x,y, label='earth Verlet')
#    x, y = np.loadtxt('euler_dt_equals_'+str(i)+'.txt', delimiter=',', unpack=True)
#    plt.plot(x,y, '--', label='earth Euler')
#    plt.xlabel('x [AU]')
#    plt.ylabel('y [AU]')
#    plt.title('Solar system, time step = {} , 1 year'.format(dts[i]))
#    plt.legend(loc='best')
#    plt.show(block=False)

#==============================================================================
# meassure time algorithms
#==============================================================================

#method_names = ['verlet', 'euler']
#system_tests = ['earth_sun','earth_sun_jupiter','all_planets','all_planets_and_moons']

#for system_test in system_tests:
#    for method_name in method_names:
#        times = np.loadtxt(method_name+'_'+system_test+'_time.txt', delimiter=',', unpack=True)
#        print('Avg time for '+ system_test +' using '+ method_name + ': ' + str(sum(times)/len(times)))


#==============================================================================
# animate solar system
#==============================================================================

#def animation_plot(method_name='verlet'):
#    fig = plt.figure()
#    list_images =[[] for i in range(0,number_bodies)]
#    list_animations = [None for i in range(0,number_bodies)]

#    for index, body in enumerate(body_list):
#        print('Loading: ', body)
#        x,y = np.loadtxt('verlet_'+body+'.txt', delimiter=',', unpack=True) 
#        print('Generating plot: ', body)
#        for i in range(0,1000,1):
#            list_images[index].append(plt.plot(x[i:i+3:1],y[i:i+3:1], label=body, color = 'blue'))
#        list_animations[index] = animation.ArtistAnimation(fig, list_images[index], interval=200)
#    plt.axis([-40, 40, -40, 40])
#    plt.xlabel('x [AU]')
#    plt.ylabel('y [AU]')
#    plt.title('Solar system')
#    plt.show()

#animation_plot()

#==============================================================================
# Study numerical error
#==============================================================================

#dts = [ 1.0 / 12, 1.0 / 52, 1.0 / 365, 1.0 / (365 * 24), 1.0 / (365 * 24 * 4) ,1.0 / (365 * 24 * 60), 1.0 / (365 * 24 * 60 * 60), 1.0 / (365 * 24 * 60 * 60 * 10), (1.0 / (365 * 24 * 60 * 60 * 10)/10)]
#x_ax = [0,1,2,3,4,5,6,7,8]
#my_xticks = ['1 month','1 week', '1 day', '1 hour', '15 minutes', '1 minute', '1 s','1/10 s', '1/100 s']
#T = 1 # max time
#day = 1.0/365

#def analytical_orbit(t):
#    x = np.cos(2*np.pi*t);
#    y = np.sin(2*np.pi*t);
#    return x, y


#def plot_error(method_name='verlet'):
#    rel_error = []
#    for index, dt in enumerate(dts):
#        x, y = np.loadtxt(method_name+'_dt_equals_'+str(index)+'.txt', delimiter=',', unpack=True)
#        N = len(x)
#        t = np.arange(0,T+dt,day)
#        x_analytical, y_analytical = analytical_orbit(t)
#        S = 0
#        for i in range(0,N):
#            r = np.array([x[i],y[i]])
#            r_analytical = np.array([x_analytical[i],y_analytical[i]])
#            S+=np.linalg.norm(r-r_analytical)
#        rel_error.append(S/N)
#        print("completed stage: ", index)
#    plt.plot(rel_error,'o-')   

#plt.figure()
#plot_error()
#plot_error('euler')
#plt.legend(['verlet','euler'])
#plt.xticks(x_ax, my_xticks)
#plt.xlabel("$dt$")
#plt.ylabel("$L^1$ relative error")
#plt.yscale('log')
#plt.title("{} years circular orbit earth-sun system".format(T))
#plt.show(block=False)
#plt.show()

#==============================================================================
# Method stability
#==============================================================================

#method_names = ['verlet', 'euler']
#for method_ in method_names:
#    t, K, U, Lx, Ly, Lz = np.loadtxt('E_L_data_'+method_+'.txt', delimiter=',', unpack=True)
#    K_new, U_new, Lx_new, Ly_new, Lz_new = [], [], [], [], []
#    tmin = min(t); tmax = max(t); N = len(t); dt = (tmax-tmin)/N
#    n = int(50 * float(N)/tmax)
#    t_new = []
#    for index in range(0,int(float(N)/n)):
#        imin = int(index*n); imax = int((index+1)*n)
#        K_new.append(max(abs((K[imin:imax]-K[0])/K[0])))
#        U_new.append(max(abs((U[imin:imax]-U[0])/U[0])))
#        Lz_new.append(max(abs((Lz[imin:imax]-Lz[0])/Lz[0])))
#        t_new.append(dt*index*n)
#    plt.figure(1)
    

#    plt.step(t_new,[K_new[0]]+K_new[:-1], label = r'Ceiling K rel error ({})'.format(method_))
#    plt.step(t_new,[U_new[0]]+U_new[:-1], linestyle='dashed', label = r'Ceiling U rel error ({})'.format(method_))

#    plt.figure(2)
#    plt.plot(t[0:-50*365],abs((Lz[0:-50*365]-Lz[0])/Lz[0]), label = 'L_z ({})'.format(method_))
#    plt.step(t_new,[Lz_new[0]]+Lz_new[:-1],label = r'Ceiling L_z ({})'.format(method_))


#plt.figure(1)
#plt.xlabel('time [yr]', fontsize=14)
##plt.ylabel(r'max$_{t}\,\,\frac{E(t)-E_0}{E_0}, \,\,\, t\in[50n,50(n+1)]$', fontsize=14)
#plt.title('Energy of Earth-Sun system', fontsize=14)
#plt.ylabel('Relative error (compared to analytical solution)', fontsize=12)
#plt.legend(loc='best')
#plt.yscale('log')


#plt.figure(2)
#plt.xlabel('time [yr]', fontsize=14)
#plt.ylabel('Relative error (compared to analytical solution)', fontsize=12)
##plt.ylabel(r'max$_{t}\,\,\frac{L_z(t)-(L_z)_0}{(L_z)_0}, \,\,\, t\in[50n,50(n+1)]$', fontsize=14)
#plt.ticklabel_format(style='sci', scilimits=(-1,3))
#plt.title('Angular momentum of Earth-Sun system', fontsize=14)
#plt.legend(loc='best')
#plt.show()

#==============================================================================
# Study escape velocity
#==============================================================================

#s_list = [ 0, 0.25, 0.5, 0.75, 1,  1.25 ]
#beta_list = [ 3, 3.25, 3.5 , 3.75, 4 ]
#for j, beta in enumerate(beta_list):
#    plt.figure(j) 
#    graph_names = []
#    for i, s in enumerate(s_list):
#        x, y = np.loadtxt('earth_escape_s_'+str(i)+'_beta_'+str(j)+'.txt', delimiter=',', unpack=True)
#        graph_names.append('$s=${}'.format(s_list[i]))
#        plt.plot(x,y)
#    plt.xlabel('x [AU]')
#    plt.ylabel('y [AU]')
#    plt.title(r"Sun-Earth system, $\beta=${}".format(beta_list[j]-1))
#    plt.legend(graph_names)
#    plt.show(block=False)
#    plt.axis([-15, 3, -5, 15])

#==============================================================================
# Earth-sun-jupiter
#==============================================================================

masses = [1,10,100,1000]




for index, mass in enumerate(masses):
    x_sun, y_sun = np.loadtxt('x_y_earth_sun_jupiter_system_sun_m_'+str(masses[index])+'.txt.', delimiter=',', unpack=True)
    plt.figure(index*100)
    x, y = np.loadtxt('x_y_earth_sun_jupiter_system_earth_m_'+str(masses[index])+'.txt.', delimiter=',', unpack=True)
    plt.plot(x-x_sun,y-y_sun)
    #plt.plot(x,y)
    x, y = np.loadtxt('x_y_earth_sun_jupiter_system_jupiter_m_'+str(masses[index])+'.txt.', delimiter=',', unpack=True)
    plt.plot(x-x_sun,y-y_sun)
    #plt.plot(x,y)
    plt.axis([-8, 8, -8, 8])
    plt.xlabel('x [AU]')
    plt.ylabel('y [AU]')
    plt.legend(['Earth','Jupiter'])
    plt.title(r"Earth-Sun-Jupiter system, $m_J$ is {} times bigger than normal".format(masses[index]))
    t, K, U, Lx, Ly, Lz = np.loadtxt('E_L_data_verlet_earth_jupiter_m_'+str(masses[index])+'.txt.', delimiter=',', unpack=True)
    E = [U[i]+K[i] for i in range(0,len(U))]
    E_new, Lx_new, Ly_new, Lz_new = [], [], [], []
    tmin = min(t); tmax = max(t); N = len(t); dt = (tmax-tmin)/N
    n = int(50 * float(N)/tmax)
    t_new = []
    for j in range(0,int(float(N)/n)):
        imin = int(j*n); imax = int((j+1)*n)
        E_new.append(max(abs((E[imin:imax]-E[0])/E[0])))
        Lz_new.append(max(abs((Lz[imin:imax]-Lz[0])/Lz[0])))
        t_new.append(dt*j*n)
    plt.figure(1)
    plt.step(t_new,[E_new[0]]+E_new[:-1], label = "Ceiling $E=K+U$ relative error {}".format(masses[index]))
    plt.figure(2)
    #plt.plot(t[0:-50*365],abs((Lz[0:-50*365]-Lz[0])/Lz[0]), label = 'L_z ({})')
    plt.step(t_new,[Lz_new[0]]+Lz_new[:-1], label = "Ceiling  $L_z$ relative error {}".format(masses[index]))

plt.figure(1)
plt.xlabel('time [yr]', fontsize=14)
plt.ylabel('relative error', fontsize=14)
plt.title('Energy of Earth-Sun-Jupiter system', fontsize=14)
plt.legend(loc='best')
plt.figure(2)
plt.yscale('log')
plt.title('Angular momentum of Earth-Sun-Jupiter system', fontsize=14)
plt.xlabel('time [yr]', fontsize=14)
plt.ylabel('relative error', fontsize=14)
plt.legend(loc='best')

#==============================================================================
# Perihelion precession
#==============================================================================


#def perihelion_precision_plot(method_name='rel', reg=True):
#    print("Loading file.")
#    if method_name == 'rel':
#        t, x, y = np.loadtxt('mercury_prehelion_precission.txt', delimiter=',', unpack=True)
#    elif method_name == 'non_rel':
#        t, x, y = np.loadtxt('mercury_prehelion_precission_non_rel.txt', delimiter=',', unpack=True)

#    print("Loading completed!")
#    r = np.sqrt(x**2+y**2)
#    rmin = np.amin(r)
#    print(rmin)
#    epsilon = 0.02
#    tolerance = rmin + epsilon
#    perihelion_indices = [i for i in range(0,len(t))]


#    # finding those r that lie within desired tolerance
#    print("Finding r within tolerance.")
#    number_items_deleted = 0
#    for i, index in enumerate(list(perihelion_indices)):
#        if r[index] > tolerance:
#            del perihelion_indices[i-number_items_deleted];
#            number_items_deleted+=1
#    print("Found all r within tolerance.")


#    #Calculating prehelion angles
#    print("Calculating prehelion angles.")
#    x_perihelion = np.array([x[i] for i in perihelion_indices])
#    y_perihelion = np.array([y[i] for i in perihelion_indices])
#    t_perihelion = np.array([t[i] for i in perihelion_indices])
#    quotient_y_by_x_perihelion = y_perihelion/x_perihelion
#    theta = np.arctan(quotient_y_by_x_perihelion)
#    theta = (3600 * 180)/np.pi*theta
#    print("Completed calculating prehelion angles.")

#    # linear regression

#    print("Finding model and plotting.")
#    slope, intercept, r_value, p_value, std_err = stats.linregress(t_perihelion,theta)
#    lin_reg_model = t_perihelion *slope + intercept
#    print(slope)

#    if method_name == 'rel':
#        plt.plot(t_perihelion, theta,'o', label = r"relativistic")
#    elif method_name == 'non_rel':
#        plt.plot(t_perihelion, theta,'o', label = r"non relativistic")
#    plt.plot(t_perihelion, lin_reg_model, label = "linear regression, slope: {}".format(slope))


#plt.figure("Perihelion precession")
#plt.title("Perihelion precession")

#perihelion_precision_plot()
#perihelion_precision_plot(method_name='non_rel')

#plt.xlabel("time [years]")
#plt.ylabel(r"perihelion angle $\theta_p$ [arc seconds]")
#plt.legend(loc='best')
#plt.show(block=False)


#=======================================================================
#Showing plots
#=======================================================================



plt.show()


#=======================================================================
# End
#=======================================================================



































#=================================================================================================
#BACKUP
#=================================================================================================

## checking that the found indices are far apart
#print("Finding smallest indices separated by one orbit")
#dt =1/(365*24*12)
#number_items_deleted = 0
#for i, index in enumerate(list(prehelion_indices[0:-1])):
#    icurr = prehelion_indices[i-number_items_deleted]
#    inext = prehelion_indices[i-number_items_deleted+1]
#    if t[inext] - t[icurr]  < dt*6000: #200 for dt=15 minutes
#        if r[icurr]>r[inext]:
#            del prehelion_indices[i-number_items_deleted]
#        else:
#            del prehelion_indices[i-number_items_deleted+1]
#        number_items_deleted+=1
#print(number_items_deleted)





#for i, index in enumerate(prehelion_indices):
#    if r[index] > tolerance:
#        #del prehelion_indices[i];
#        prehelion_indices.remove(prehelion_indices[i])









# checking that the found indices are far apart
#indices_far_apart = False
#while indices_far_apart == False:
#    N = len(prehelion_indices)
#    for i in range(0,N-1):
#        index_i = prehelion_indices[i]
#        index_i_plus_1 = prehelion_indices[i+1]
#        if i == N-2:
#            indices_far_apart = True;
#        if index_i_plus_1-index_i  < 200:
#            if r[index_i]>r[index_i_plus_1]:
#                del prehelion_indices[i]
#                break
#            else:
#                del prehelion_indices[i+1]
#                break
#        if i == N-2:
#            indices_far_apart = True;











## finding r (by index) such that r is a local minimum
#print("Finding local minima.")
#for index in range(1,len(r)-1):
#    rprev = r[index-1]
#    rcurr = r[index]
#    rnext = r[index+1]
#    if rcurr <  rprev and rcurr < rnext:
#        prehelion_indices.append(index)
#print("Finding local minima completed!")