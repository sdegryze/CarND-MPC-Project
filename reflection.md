# The Model

The Model Predictive Control approach implemented is using a kinematic model. The kinematic model uses 4 state variables:

 * vehicle x position: $x$
 * vehicle y position: $y$
 * vehicle speed: $v$
 * vehicle orientation (psi): $\psi$

The state is being impacted by 2 actuators (control variables):

 * throttle/break or acceleration: $a$
 * car steering angle: $\delta$

The kinematic model describes how the state of the vehicle at time $t+1$ can be calculated based on the state of the vehicle at time $t$ using the following equations:

$\begin{align}
x_{t+1} & = x_t + v_t * cos(\psi_t) * dt \nonumber \\
y_{t+1} & = y_t + v_t * sin(\psi_t) * dt \nonumber \\
\psi_{t+1} & = \psi_t + \frac{v_t}{L_f} * \delta_t * dt \nonumber \\
v_{t+1} & = v_t + a_t * dt \nonumber \\
\end{align}$

In addition, the cross-track error and orientation error are calculated as follows:

$\begin{align}
cte_{t+1} & = f(x_t) - y_t + v_t * sin(e\psi_t) * dt \nonumber \\
e\psi_{t+1} &= \psi_t - {\psi}des_t + \frac{v_t}{L_f} * \delta_t * dt \nonumber
\end{align}$

The cost was calculated as following:

$\begin{align}
J = \sum_{t=1}^{N}cte_t^2 + e\psi_t^2 + (v_t - v_{ref})^2 + 35000\delta_t^2 + 5000(\delta_t - \delta_{t-1})^2 + 10(a_t - a_{t-1})^2
\end{align}$

The main value I tuned was the factor for $\delta_t^2$, i.e., 35000. High values made the car slow in reacting to curves. Lower values made the car take too sharp turns. I found that 35000 was a good middle ground where the vehicle was able to drive fast while staying on the driveable part of the road.

# Timestep Length and Elapsed Duration (N & dt)

I ended up with an N of 30 and a duration of a timestep dt of 0.025 sec. This gave me a look-ahead horizon of 0.75 sec. This was also about the duration I was passing the waypoints that are provided by the simulator to the model each time. If I increased the N so that the horizon was much longer than 0.75 sec, the polynomial that was being fitted would extend far beyond the waypoints provided and often not follow the curve of the road, leading to bias. I started off with a timestep of 0.1 sec, but found that this was too coarse to follow the track with enough "reactivity", especially at higher speeds (75 mph). Reducing the timestep to 0.025 is perhaps a little overkill, but definitely made a noticeable improvement.

# Polynomial Fitting and MPC Preprocessing

Before a 3rd degree polynomial fit was executed, I transformed the waypoint coordinates from the map coordinate system (as provided by the simulator) into the car coordinate system. No other pre-processing was done.

# Model Predictive Control with Latency

I first handled the latency with a higher value of the factor for $\delta_t^2$ in the cost calculation (more towards 80000). That worked fairly well. I also tried to select not the first estimate of the speed and throttle in the MPC solution, but rather the 4th (which is 100 ms into the future). That did not work well. The car became very sluggish and not reactive at all. Finally, I explicitly incorporated a 100 ms latency by assuming the vehicle was 100 ms further than the values I received from the simulator. Subsequently, all calculations would be for this new position, so that commands corresponding for this position reach the actuators at the time they were simulated for.