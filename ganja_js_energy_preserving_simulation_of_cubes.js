/** May The Forque Be With You **/

/* A dimension independent approach to dynamics */
/* Leo Dorst & Steven De Keninck.               */
/* 20220102                                     */

/* Adapted by vlad2305m in 2025                 */

// Set the number of dimensions. 
window.d = 3;

// Create a d-dimensional geometric algebra over the reals.
Algebra(d, 0, 1, ()=>{
  
  // Create a hypercube (square in 2D, cube in 3D, ...)
  // Start by defining its vertices.
  var p = [...Array(2**d)]                         
          .map((x,i)=>i.toString(2))              // [0, 1, 10, 11, 100, ...]
          .map(x=>('00000'+x).slice(-d))          // [000, 001, 010, 011, ...]
          .map(x=>x.split('').map(x=>x-0.5))      // [[-0.5, -0.5, -0.5], [-0.5, -0.5, 0.5], ...]
          .map(x=>!(1e0 + x*[1e1,1e2,1e3,1e4,1e5]));   // PGA points are dual vectors.
  
  // Now consider all vertex pairs and create edges for those
  // pairs that differ only in one coordinate.
  var e = p.map((a,i)=>p.map((b,j)=>
            i<=j||(i^j)&(i^j-1)?0:[a,b]   // note that &,^ here are bitwise ops since i,j are integer
          )).flat();   
          
  //  Attachement points for Hooke's law: 1 in body frame and 1 in world frame
  var spring1 = [p[2**d-1], 1*p[2**d-1]];              // last point of hypercube.
  var spring2 = [p[2**(d-1)], 1*p[2**(d-1)]];
  
  var floor  = (1e2) + 3e0;           // tilted hyperplane
  
  // Physics state.    
  var states = [[(Math.E**(.1e12)*(1+0.5*1e02+0.5e03)), 0.1e01+0.1e12, spring1,[1]],
                [(Math.E**(.1e12)*(1-0.5*1e02-0.5e03)), 0.1e01+0.1e12, spring2,[0]]];
  
  var g = -9.81e02;
  var k = 16;
  var repK = 100;
  var Hooke = (M, a) => k*( (~M>>>a[1]) & a[0] );
  // Our forces (expressed in the body frame, as forque lines)
  var F = ([M,B], s)=>{
        var Gravity = !(~M >>> g);
        var Damping = -.5*!B;
        return Gravity + Damping*0 + Hooke(M, s);
      };
      
  var Frep = (M, others) => {
    let forque = 0e0;
    for (const other of others) {
      const [M1] = states[other];
      let F1 = -((~M >>> (M1 >>> !1e0)) & !1e0);
      F1 = repK * F1 / (F1.Length**3);
      forque = forque + F1;
    }
    return forque;
  }
  var Erep = (M, others) => {
    let E = 0;
    for (const other of others) {
      const [M1] = states[other];
      let d1 = ((~M >>> (M1 >>> !1e0)) & !1e0).Length;
      let E1 = repK / d1;
      E = E + E1;
    }
    return E;
  }
  
  var stepQ = ([M,B,...extra], dt) => {
    var diffq = dt * -0.5*M*B;
    M = M + diffq;
    return [M.Normalized,B.Grade(2),...extra];
  }
  const dims = [1e01, 1e02, 1e03, 1e04, 1e05, 1e06, 1e07];
  var stepQrot = ([M,B,...extra], dt, iterdir) => {
    // UnDual of the commutator only gives us components in e0i, which is convenient for decomposition
    for (let i = 0; i < d; i++) {
      var j = iterdir ? i : d-i-1;
      var comm = (-0.5*(B.Dual*B-B*B.Dual)).UnDual;
      var diffp = dt * (comm[j+2+d]*dims[j]); // project to dims[i]
      //if (c%100===0) console.log(""+(i)+" "+diffp/dt) // check formula is correct
      B = B + diffp;
    }
    return [M.Normalized,B.Grade(2),...extra];
  }
  var stepP = ([M,B,s,others], dt) => {
    let forque = F([M,B], s) + Frep(M, others);
    var diffp = dt * forque.UnDual;
    B = B + diffp;
    return [M,B,s,others];
  }
  
  var vec = (fun, states, ...args) => {
    const states1=[]
    for (const state of states) 
      states1.push(fun(state, ...args));
    return states1;
  }

  var c = 0;
  var Ebounds = [1/0, -1/0];
  var TEbounds = [0,0];
  // Render
  return this.graph(()=>{
    c++;
    var T = performance.now()/1000;
    const dt = 1/600;
    for (var i=0; i<10; ++i) {
      states = vec(stepQ, states, dt/2);
      states = vec(stepP, states, dt);
      states = vec(stepQ, states, dt/2);
      states = vec(stepQrot, states, dt, i%2);
    }
    var E = 0;
    for (const [M,B,s,others] of states) {
      E = E + ((B&!B) + (Hooke(M, s).Length**2)/k - 2*(1-2*(d%2))*(((!1e02|(M>>>!1e0))*!g)&!1e0)*(!g).Normalized);
      E = E + Erep(M, others);
    }
    if (E > Ebounds[1]) {
      Ebounds[1] = E;
      TEbounds[1] = T;
    }
    if (E < Ebounds[0]) {
      Ebounds[0] = E;
      TEbounds[0] = T;
    }
    return [
      "Emax&hairsp;="+Ebounds[1] + ` ${Math.floor((T-TEbounds[1])/60)}min ${Math.floor(T-TEbounds[1])%60}s`,
      "E&emsp;&emsp;&thinsp;="+E,
      "Emin ="+Ebounds[0] + ` ${Math.floor((T-TEbounds[0])/60)}min ${Math.floor(T-TEbounds[0])%60}s`,
      0x009977, floor*2,
      ...states.map( ([M,B,s]) => [
        0xccccff, ...p.map(x=>M >>> [x,x-0.05*(x*B-B*x)]),
        0,        ...M >>> e,
        0x007799, s[1], [s[1], M>>>s[0]]
      ]).flat(1),
      ]
  
  },{lineWidth:3,pointRadius:2,animate:1,scale:0.75});
})
