#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double velocity = j[1]["speed"];
     	    double delta = j[1]["steering_angle"];
		      double accel = j[1]["throttle"];
          

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
	 	  // Vectors for storing x/y coords for the waypoints
          vector<double> wayPt_x;
          vector<double> wayPt_Y;


          for (int i = 0; i < ptsx.size(); i++) {
          
       // Push calculated values back to Waypoint X/Y Vectors. Homogeneous TM used
       // https://www.miniphysics.com/coordinate-transformation-under-rotation.html
		   // Transforms the perspective for the waypoint to be from the vehicles point of view.
          double x_Shift = ptsx[i] - px;
          double y_Shift = ptsy[i] - py;
          wayPt_x.push_back( x_Shift * cos(-psi) - y_Shift * sin(-psi) );
          wayPt_Y.push_back( x_Shift * sin(-psi) + y_Shift * cos(-psi) );
          }

          // Set steeting and Throttle values
          double steer_value = j[1]["steering_angle"];
          double throttle_value = j[1]["throttle"];

          //cout << 'The Steering Angle is' << steer_value;
          //cout << 'The Throttle input is' << throttle_value;

          // In order to type cast vector<doubles> to 
          double* pointerX = &wayPt_x[0];
          double* pointerY = &wayPt_Y[0];

          Eigen::Map<Eigen::VectorXd> wayPt_x_Eig( pointerX, 6 );        
          Eigen::Map<Eigen::VectorXd> wayPt_y_Eig( pointerY, 6 );

          //cout << sdata << endl;
          //Lesson 20 Fitting Polynomials to the waypoints
          // From Lesson 19-8 usually a 3rd degree poly is used fits most roads.
          auto coeffs = polyfit( wayPt_x_Eig, wayPt_y_Eig, 3 );
          // Cross track and EPSI error calculation.
          //                    ^Cars Point^ - ^Polynomial^
          //double cte = polyeval(coeffs, x) - y;

          // This provides an estimated CTE.  THe zero is due to the 
          // Shifting of the reference frame on line 113 and 114.  This
          // gives us a y value.
          double cte = polyeval( coeffs, 0 );
          // Simplified due to the px = 0, py = 0 and psi = 0 assumption.
          // -atan( a + bx + cx^2 )------> -atan( a ) ******  [bc x = 0] 
          double epsi = -atan( coeffs[1] );  

		 		  
          const double del_t = 0.1;
          const double Lf = 2.67;
		      const double px_t = 0.0 + velocity * del_t;
		      const double py_t = 0.0;
		      const double psi_t = 0.0 - ( velocity * delta * del_t )/ Lf;
		      const double vel_t = velocity + accel * del_t;
		      const double cte_t = cte + velocity * sin(epsi) * del_t;
    		  const double epsi_t = epsi - ( velocity * del_t * delta ) / Lf;
		   
          // 
          Eigen::VectorXd state(6);
          state << px_t, py_t, psi_t, vel_t, cte_t, epsi_t;
          auto vars = mpc.Solve(state, coeffs);
          steer_value = vars[0];
          throttle_value = vars[1];

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value/(deg2rad(25));
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          for (int i = 2; i < vars.size(); i ++) {
            if (i%2 == 0) {
              mpc_x_vals.push_back(vars[i]);
            }
            else {
              mpc_y_vals.push_back(vars[i]);
            }
          }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          for (double i = 0; i < wayPt_x.size(); i ++){
            next_x_vals.push_back(wayPt_x[i]);
            next_y_vals.push_back(wayPt_Y[i]);
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}