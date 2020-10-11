// Scaling variables
var scaleX;
var scaleY;
var sW1;
var sW2;

// Ellipse parameters (for testing) 
var pVal;
var qVal;
var wVal;
var hVal;
var tVal;
var movableEllipse;

// Initial iterate (for testing)
var xInitial;
var yInitial;

// Number of points in orbit (for testing)
var nVal;

// Orbits in global space
var ponceletOrbit = [];
var ponceletInnerOrbit = [];

// Durations in global space
var angularDurations = [];
var euclideanDurations = [];
var innerAngularDurations = [];
var innerEuclideanDurations = [];

// Use this to export orbit data 
var dataExport = false;

// Points of intersection on the circumscribing square
var squareIntersections = [];

var testPoints = [];

var unitCircleOn;
var ellipseOn;
var squareOn;
var ponceletOrbitOn;
var linesOn;
var curvesOn;
var beziersOn;
var bezier0On, bezier1On, bezier2On, bezier3On, bezier4On, bezier5On;

export function setup(audioData, sendData) {
  createCanvas(600, 600);

  scaleX = 200;
  scaleY = 200;
  sW1 = 2/scaleX;
  sW2 = 0.5/scaleX;

  // pVal = 0.4;
  // qVal = 0.2;
  // wVal = 0.5;
  // hVal = 0.25;
  // tVal = PI/6.0;

  // pentagon p=q=t=0, w=h=cos(PI/5)
  pVal = 0.35;
  qVal = 0.45;
  wVal = 0.2; //cos(PI/5.0); // 0.2;
  hVal = 0.2; //cos(PI/5.0); // 0.2;
  tVal = 0.0; 

  movableEllipse = new MovableEllipse(pVal,qVal,wVal,hVal,tVal);

  xInitial = cos(3*PI/6.0);
  yInitial = sin(3*PI/6.0);

  nVal = 150;

  unitCircleOn = true;
  ellipseOn = true;
  squareOn = false;
  ponceletOrbitOn = false;
  linesOn = false;
  curvesOn = false;
  beziersOn = true;
  bezier0On = false;
  bezier1On = false;
  bezier2On = true;
  bezier3On = false;
  bezier4On = false; 
  bezier5On = false; 

  ponceletOrbit = computeOrbit(xInitial, yInitial, movableEllipse.getp(), movableEllipse.getq(), movableEllipse.getw(), movableEllipse.geth(), movableEllipse.gett(), nVal);
  ponceletInnerOrbit = computeInnerOrbit(ponceletOrbit, movableEllipse.getp(), movableEllipse.getq(), movableEllipse.getw(), movableEllipse.geth(), movableEllipse.gett(), nVal);

  angularDurations = computeAngularDurations(ponceletOrbit);
  euclideanDurations = computeEuclideanDurations(ponceletOrbit);
  innerAngularDurations = computeAngularDurations(ponceletInnerOrbit);
  innerEuclideanDurations = computeEuclideanDurations(ponceletInnerOrbit);


  if (dataExport) {
  	saveJSON(ponceletOrbit, 'orbitData.json');
  	saveJSON(angularDurations, 'angularDurations.json');
  	saveJSON(euclideanDurations, 'euclideanDurations.json');
  }

  squareIntersections = computeSquareIntersections(ponceletOrbit);
	sendData({
		pVal,
		qVal,
		ponceletOrbit,
		ponceletInnerOrbit,
		angularDurations ,
		euclideanDurations,
		innerAngularDurations,
		innerEuclideanDurations,
	})

	movableEllipse.onDrag(sendData)
}

export function draw(audioData) {
	background(255);
  
	scale(scaleX,-scaleY);
	translate((width/2)/scaleX,-(height/2)/scaleY);

	if (squareOn) {
		drawCicumscribingSquare();
	}

	if (unitCircleOn) {
		drawUnitCircle();
	} 

	if (ellipseOn) {
		movableEllipse.display();
		movableEllipse.drag(mouseX, mouseY);
	}

	if (ponceletOrbitOn) {
		drawOrbit(ponceletOrbit);
		// drawInnerOrbit(ponceletInnerOrbit);
	}

	if (linesOn) {
		drawSquareIntersectionsLine(squareIntersections);
	}

	if (curvesOn) {		
		drawSquareIntersectionsCurve(squareIntersections);
	}

	if (beziersOn) {
		if (bezier0On) {		
			drawSquareIntersectionsBezier(squareIntersections);
		}
		if (bezier1On) {
			drawBezierFlower1(ponceletOrbit,ponceletInnerOrbit);
		}
		if (bezier2On) {
			drawBezierFlower2(ponceletOrbit,ponceletInnerOrbit);
		}
		if (bezier3On) {
			drawBezierFlower3(ponceletOrbit,ponceletInnerOrbit);
		}
		if (bezier4On) {
			drawBezierFlower4(ponceletOrbit,ponceletInnerOrbit);
		}
		if (bezier5On) {
			drawBezierFlower5(ponceletOrbit,ponceletInnerOrbit);
		}		
	}

}

class MovableEllipse {	

  constructor(p, q, w, h, t) {
    this.p = p;
    this.q = q;
    this.h = h;
    this.w = w;
    this.h = h;
    this.t = t;
    this.fillColour = color(255,255,255);
    this.fillOn = false;
    this.dragging = false;
    this.offsetp = 0;
    this.offsetq = 0;
    this.transformedMouseXY = [];
    this.pFixed = 0;
    this.qFixed = 0;
    this.tFixed = 0;
    this.pInside = 0;
    this.qInside = 0;
    this.alpha = 0;
  }

  clicked(mx, my) {
  	this.transformedMouseXY = transformCoordinatesEllipse(mx, my, this.p, this.q, this.t);
  	if (this.transformedMouseXY[0] > -this.w && this.transformedMouseXY[0] < this.w && this.transformedMouseXY[1] > -this.h && this.transformedMouseXY[1] < this.h) {
		this.dragging = true;
		this.offsetp = this.p - this.transformedMouseXY[0];
		this.offsetq = this.q - this.transformedMouseXY[1];
		this.pFixed = this.p;
		this.qFixed = this.q;
		this.tFixed = this.t;
  	}
  }

  drag(mx, my) {
		let [pX, pY] = this._prevMouseXY || [];
		if(pX === mx && pY === my) return;
  	this.transformedMouseXY = transformCoordinatesEllipse(mx, my, this.pFixed, this.qFixed, this.tFixed);
  	if (ellipseInsideCircle(this.p, this.q, this.w, this.h)) { 	
	  	if (this.dragging) {
	    	this.p = this.transformedMouseXY[0] + this.offsetp;
	    	this.q = this.transformedMouseXY[1] + this.offsetq;
	    	this.pInside = this.p;
	    	this.qInside = this.q;
	    	ponceletOrbit = computeOrbit(xInitial, yInitial, this.p, this.q, this.w, this.h, this.t, nVal);
	    	ponceletInnerOrbit = computeInnerOrbit(ponceletOrbit, this.p, this.q, this.w, this.h, this.t, nVal);
	    	angularDurations = computeAngularDurations(ponceletOrbit);
			euclideanDurations = computeEuclideanDurations(ponceletInnerOrbit);
	    	innerAngularDurations = computeAngularDurations(ponceletInnerOrbit);
			innerEuclideanDurations = computeEuclideanDurations(ponceletOrbit);			
	    	squareIntersections = computeSquareIntersections(ponceletOrbit);
				if(frameCount % 10 === 0) {
					this._onDrag({
						pVal,
						qVal,
						ponceletOrbit,
						ponceletInnerOrbit,
						angularDurations ,
						euclideanDurations,
						innerAngularDurations,
						innerEuclideanDurations,
					})
				}
		}
	} else {
	    this.alpha = atan2(this.qInside,this.pInside);
		this.p = (1-(max(this.w,this.h)+0.0001))*cos(this.alpha);
		this.q = (1-(max(this.w,this.h)+0.0001))*sin(this.alpha);	
		this.t = this.tFixed;
		this.dragging = false;
	}
		this._prevMouseXY = [mx, my];
  }

  stopDragging() {
  	this.dragging = false;
  }

	onDrag(cb) {
		this._onDrag = cb
	}

  getp() {
  	return this.p;
  }

  getq() {
  	return this.q;
  }

  getw() {
  	return this.w;
  }

  geth() {
  	return this.h;
  }

  gett() {
  	return this.t;
  } 

  setFill(r, g, b) {
  	this.fillOn = true;
  	this.fillColour = color(r,g,b);
  }

  clearFill() {
  	this.fillOn = false;
  }     

  display() {
	push();
		stroke(0,0,0);
		if (this.fillOn) {
			fill(this.fillColour);
		} else {
			noFill();
		}
		translate(this.p,this.q);
		rotate(this.t);
		strokeWeight(sW2);
		ellipse(0,0,2*this.w,2*this.h);
	pop();	  	
  }
}

// Given a point (x,y) on the the unit circle, generate it's 
// Poncelet iterate (x',y') on the unit circle. 
// (p,q) is the centre of the ellipse
// w and h are the semi major and semi minor axes of the ellipse 
// respectively
// t is the angle of rotation of the ellipse (in radians) 
function ponceletIterate(x, y, p, q, w, h, t) {
	let iterate = [];
	let dummyIterate = [0,1];
	let sqrtArg = -(pow(w,4)*pow(h,4)*(pow(w,2)*pow(h,2)-(pow(q,2)*pow(w,2)+pow(p,2)*pow(h,2)-2*p*pow(h,2)*x+pow(h,2)*pow(x,2)-2*q*pow(w,2)*y+pow(w,2)*pow(y,2))*pow(cos(t),2)-(pow(p,2)*pow(w,2)+pow(q,2)*pow(h,2)-2*p*pow(w,2)*x+pow(w,2)*pow(x,2)-2*q*pow(h,2)*y+pow(h,2)*pow(y,2))*pow(sin(t),2)+p*q*pow(w,2)*sin(2*t)-p*q*pow(h,2)*sin(2*t)-q*pow(w,2)*x*sin(2*t)+q*pow(h,2)*x*sin(2*t)-p*pow(w,2)*y*sin(2*t)+p*pow(h,2)*y*sin(2*t)+pow(w,2)*x*y*sin(2*t)-pow(h,2)*x*y*sin(2*t)));
	if (sqrtArg>=0) {
		iterate[0] = (-(pow(p,4)*pow(w,2)*pow(h,2)*x)+pow(q,4)*pow(w,2)*pow(h,2)*x+pow(p,2)*pow(w,4)*pow(h,2)*x-pow(q,2)*pow(w,4)*pow(h,2)*x+pow(p,2)*pow(w,2)*pow(h,4)*x-pow(q,2)*pow(w,2)*pow(h,4)*x+4*pow(p,3)*pow(w,2)*pow(h,2)*pow(x,2)-2*p*pow(w,4)*pow(h,2)*pow(x,2)-2*p*pow(w,2)*pow(h,4)*pow(x,2)-6*pow(p,2)*pow(w,2)*pow(h,2)*pow(x,3)+pow(w,4)*pow(h,2)*pow(x,3)+pow(w,2)*pow(h,4)*pow(x,3)+4*p*pow(w,2)*pow(h,2)*pow(x,4)-pow(w,2)*pow(h,2)*pow(x,5)-2*pow(p,3)*q*pow(w,2)*pow(h,2)*y-2*p*pow(q,3)*pow(w,2)*pow(h,2)*y+2*p*q*pow(w,4)*pow(h,2)*y+2*p*q*pow(w,2)*pow(h,4)*y+6*pow(p,2)*q*pow(w,2)*pow(h,2)*x*y-2*pow(q,3)*pow(w,2)*pow(h,2)*x*y-6*p*q*pow(w,2)*pow(h,2)*pow(x,2)*y+2*q*pow(w,2)*pow(h,2)*pow(x,3)*y+2*pow(p,3)*pow(w,2)*pow(h,2)*pow(y,2)+6*p*pow(q,2)*pow(w,2)*pow(h,2)*pow(y,2)-2*p*pow(w,4)*pow(h,2)*pow(y,2)-2*p*pow(w,2)*pow(h,4)*pow(y,2)-6*pow(p,2)*pow(w,2)*pow(h,2)*x*pow(y,2)+pow(w,4)*pow(h,2)*x*pow(y,2)+pow(w,2)*pow(h,4)*x*pow(y,2)+6*p*pow(w,2)*pow(h,2)*pow(x,2)*pow(y,2)-2*pow(w,2)*pow(h,2)*pow(x,3)*pow(y,2)-6*p*q*pow(w,2)*pow(h,2)*pow(y,3)+2*q*pow(w,2)*pow(h,2)*x*pow(y,3)+2*p*pow(w,2)*pow(h,2)*pow(y,4)-pow(w,2)*pow(h,2)*x*pow(y,4)-4*p*q*x*sqrt(sqrtArg)+4*q*pow(x,2)*sqrt(sqrtArg)+2*pow(p,2)*y*sqrt(sqrtArg)-2*pow(q,2)*y*sqrt(sqrtArg)-2*pow(x,2)*y*sqrt(sqrtArg)+4*q*pow(y,2)*sqrt(sqrtArg)-2*pow(y,3)*sqrt(sqrtArg)+(pow(w,2)-pow(h,2))*sin(2*t)*(pow(p,2)*pow(w,2)*pow(h,2)*y+pow(q,2)*pow(w,2)*pow(h,2)*y-pow(w,4)*pow(h,2)*y-pow(w,2)*pow(h,4)*y-2*p*pow(w,2)*pow(h,2)*x*y+pow(w,2)*pow(h,2)*pow(x,2)*y-2*q*pow(w,2)*pow(h,2)*pow(y,2)+pow(w,2)*pow(h,2)*pow(y,3)+2*x*sqrt(sqrtArg))+(pow(w,2)-pow(h,2))*cos(2*t)*(pow(p,2)*pow(w,2)*pow(h,2)*x+pow(q,2)*pow(w,2)*pow(h,2)*x-pow(w,4)*pow(h,2)*x-pow(w,2)*pow(h,4)*x-2*p*pow(w,2)*pow(h,2)*pow(x,2)+pow(w,2)*pow(h,2)*pow(x,3)-2*q*pow(w,2)*pow(h,2)*x*y+pow(w,2)*pow(h,2)*x*pow(y,2)-2*y*sqrt(sqrtArg)))/(pow(w,2)*pow(h,2)*(pow(p,4)+2*pow(p,2)*pow(q,2)+pow(q,4)+pow(w,4)-2*pow(w,2)*pow(h,2)+pow(h,4)-4*pow(p,3)*x-4*p*pow(q,2)*x+6*pow(p,2)*pow(x,2)+2*pow(q,2)*pow(x,2)-4*p*pow(x,3)+pow(x,4)-4*pow(p,2)*q*y-4*pow(q,3)*y+8*p*q*x*y-4*q*pow(x,2)*y+2*pow(p,2)*pow(y,2)+6*pow(q,2)*pow(y,2)-4*p*x*pow(y,2)+2*pow(x,2)*pow(y,2)-4*q*pow(y,3)+pow(y,4)-2*(pow(w,2)-pow(h,2))*(pow(p,2)-pow(q,2)-2*p*x+pow(x,2)+2*q*y-pow(y,2))*cos(2*t)-4*(pow(w,2)-pow(h,2))*(p-x)*(q-y)*sin(2*t)));
		iterate[1] = (-2*pow(p,3)*q*pow(w,2)*pow(h,2)*x-2*p*pow(q,3)*pow(w,2)*pow(h,2)*x+2*p*q*pow(w,4)*pow(h,2)*x+2*p*q*pow(w,2)*pow(h,4)*x+6*pow(p,2)*q*pow(w,2)*pow(h,2)*pow(x,2)+2*pow(q,3)*pow(w,2)*pow(h,2)*pow(x,2)-2*q*pow(w,4)*pow(h,2)*pow(x,2)-2*q*pow(w,2)*pow(h,4)*pow(x,2)-6*p*q*pow(w,2)*pow(h,2)*pow(x,3)+2*q*pow(w,2)*pow(h,2)*pow(x,4)+pow(p,4)*pow(w,2)*pow(h,2)*y-pow(q,4)*pow(w,2)*pow(h,2)*y-pow(p,2)*pow(w,4)*pow(h,2)*y+pow(q,2)*pow(w,4)*pow(h,2)*y-pow(p,2)*pow(w,2)*pow(h,4)*y+pow(q,2)*pow(w,2)*pow(h,4)*y-2*pow(p,3)*pow(w,2)*pow(h,2)*x*y+6*p*pow(q,2)*pow(w,2)*pow(h,2)*x*y-6*pow(q,2)*pow(w,2)*pow(h,2)*pow(x,2)*y+pow(w,4)*pow(h,2)*pow(x,2)*y+pow(w,2)*pow(h,4)*pow(x,2)*y+2*p*pow(w,2)*pow(h,2)*pow(x,3)*y-pow(w,2)*pow(h,2)*pow(x,4)*y+4*pow(q,3)*pow(w,2)*pow(h,2)*pow(y,2)-2*q*pow(w,4)*pow(h,2)*pow(y,2)-2*q*pow(w,2)*pow(h,4)*pow(y,2)-6*p*q*pow(w,2)*pow(h,2)*x*pow(y,2)+6*q*pow(w,2)*pow(h,2)*pow(x,2)*pow(y,2)-6*pow(q,2)*pow(w,2)*pow(h,2)*pow(y,3)+pow(w,4)*pow(h,2)*pow(y,3)+pow(w,2)*pow(h,4)*pow(y,3)+2*p*pow(w,2)*pow(h,2)*x*pow(y,3)-2*pow(w,2)*pow(h,2)*pow(x,2)*pow(y,3)+4*q*pow(w,2)*pow(h,2)*pow(y,4)-pow(w,2)*pow(h,2)*pow(y,5)+2*pow(p,2)*x*sqrt(sqrtArg)-2*pow(q,2)*x*sqrt(sqrtArg)-4*p*pow(x,2)*sqrt(sqrtArg)+2*pow(x,3)*sqrt(sqrtArg)+4*p*q*y*sqrt(sqrtArg)-4*p*pow(y,2)*sqrt(sqrtArg)+2*x*pow(y,2)*sqrt(sqrtArg)-(pow(w,2)-pow(h,2))*cos(2*t)*(pow(p,2)*pow(w,2)*pow(h,2)*y+pow(q,2)*pow(w,2)*pow(h,2)*y-pow(w,4)*pow(h,2)*y-pow(w,2)*pow(h,4)*y-2*p*pow(w,2)*pow(h,2)*x*y+pow(w,2)*pow(h,2)*pow(x,2)*y-2*q*pow(w,2)*pow(h,2)*pow(y,2)+pow(w,2)*pow(h,2)*pow(y,3)+2*x*sqrt(sqrtArg))+(pow(w,2)-pow(h,2))*sin(2*t)*(pow(p,2)*pow(w,2)*pow(h,2)*x+pow(q,2)*pow(w,2)*pow(h,2)*x-pow(w,4)*pow(h,2)*x-pow(w,2)*pow(h,4)*x-2*p*pow(w,2)*pow(h,2)*pow(x,2)+pow(w,2)*pow(h,2)*pow(x,3)-2*q*pow(w,2)*pow(h,2)*x*y+pow(w,2)*pow(h,2)*x*pow(y,2)-2*y*sqrt(sqrtArg)))/(pow(w,2)*pow(h,2)*(pow(p,4)+2*pow(p,2)*pow(q,2)+pow(q,4)+pow(w,4)-2*pow(w,2)*pow(h,2)+pow(h,4)-4*pow(p,3)*x-4*p*pow(q,2)*x+6*pow(p,2)*pow(x,2)+2*pow(q,2)*pow(x,2)-4*p*pow(x,3)+pow(x,4)-4*pow(p,2)*q*y-4*pow(q,3)*y+8*p*q*x*y-4*q*pow(x,2)*y+2*pow(p,2)*pow(y,2)+6*pow(q,2)*pow(y,2)-4*p*x*pow(y,2)+2*pow(x,2)*pow(y,2)-4*q*pow(y,3)+pow(y,4)-2*(pow(w,2)-pow(h,2))*(pow(p,2)-pow(q,2)-2*p*x+pow(x,2)+2*q*y-pow(y,2))*cos(2*t)-4*(pow(w,2)-pow(h,2))*(p-x)*(q-y)*sin(2*t)));		
	}

	if (iterate.length>0) {
		return iterate;
	} else {
		return dummyIterate;
	}
	
}

function ponceletInnerIterate(x0, y0, x1, y1, p, q, w, h, t) {
	let iterate = [];
	iterate[0]=((p*pow(h,2)*pow(x0-x1,2)+pow(w,2)*(y0-y1)*(q*(x0-x1)+x1*y0-x0*y1))*pow(cos(t),2)-(q*(-(pow(h,2)*pow(x0-x1,2))+pow(w,2)*(pow(x0,2)+pow(x1,2)))+(pow(w,2)-pow(h,2))*(x0-x1)*(x1*y0+p*(y0-y1)-x0*y1))*cos(t)*sin(t)+sin(t)*(2*q*pow(w,2)*x0*x1*cos(t)+(p*pow(w,2)*pow(x0-x1,2)+pow(h,2)*(y0-y1)*(q*(x0-x1)+x1*y0-x0*y1))*sin(t)))/((pow(h,2)*pow(x0-x1,2)+pow(w,2)*pow(y0-y1,2))*pow(cos(t),2)+sin(t)*(-2*(pow(h,2)*(x1*y0+x0*y1)+pow(w,2)*(x0*y0+x1*y1))*cos(t)+(pow(w,2)*pow(x0-x1,2)+pow(h,2)*pow(y0-y1,2))*sin(t))+(pow(w,2)*(x1*y0+x0*y1)+pow(h,2)*(x0*y0+x1*y1))*sin(2*t));
	iterate[1]=((p*pow(h,2)*(x0-x1)*(y0-y1)+q*pow(w,2)*pow(y0-y1,2)+pow(h,2)*(x0-x1)*(-(x1*y0)+x0*y1))*pow(cos(t),2)-(q*(pow(w,2)-pow(h,2))*(x0-x1)*(y0-y1)-p*pow(h,2)*pow(y0-y1,2)-(pow(w,2)-pow(h,2))*(y0-y1)*(x1*y0-x0*y1)+p*pow(w,2)*(pow(y0,2)+pow(y1,2)))*cos(t)*sin(t)+sin(t)*(2*p*pow(w,2)*y0*y1*cos(t)+(p*pow(w,2)*(x0-x1)*(y0-y1)+q*pow(h,2)*pow(y0-y1,2)+pow(w,2)*(x0-x1)*(-(x1*y0)+x0*y1))*sin(t)))/((pow(h,2)*pow(x0-x1,2)+pow(w,2)*pow(y0-y1,2))*pow(cos(t),2)+sin(t)*(-2*(pow(h,2)*(x1*y0+x0*y1)+pow(w,2)*(x0*y0+x1*y1))*cos(t)+(pow(w,2)*pow(x0-x1,2)+pow(h,2)*pow(y0-y1,2))*sin(t))+(pow(w,2)*(x1*y0+x0*y1)+pow(h,2)*(x0*y0+x1*y1))*sin(2*t));

	return iterate;
}

function drawUnitCircle() {
	push();
		stroke(0,0,0);
		noFill();
		strokeWeight(sW2);
		circle(0, 0, 2);
	pop();
}

// See comment above definition of ponceletIterate function
// for explanation of p, q, w, h, t
function drawEllipse(p,q,w,h,t) {
	push();
		stroke(255,0,0);
		noFill();
		translate(p,q);
		rotate(t);
		strokeWeight(sW1);
		ellipse(0,0,2*w,2*h);
	pop();	
}

// See comment above definition of ponceletIterate function
// for explanation of p, q, w, h, t
// n is the number of points in the orbit
function computeOrbit(x, y, p, q, w, h, t, n) {
	let orbit = [];
	let iterate = [];
	let iterateOld = [];
	let iterateNew = [];
	iterateOld[0] = x;
	iterateOld[1] = y;
	orbit[0] = [];
	orbit[0][0] = iterateOld[0];
	orbit[0][1] = iterateOld[1];
	for (let i=1; i<n; i++) {
		orbit[i] = [];
		iterate = ponceletIterate(iterateOld[0], iterateOld[1], p, q, w, h, t);
		iterateNew[0] = iterate[0];
		iterateNew[1] = iterate[1];
		orbit[i][0] = iterateNew[0];
		orbit[i][1] = iterateNew[1];
		iterateOld[0] = iterateNew[0];
		iterateOld[1] = iterateNew[1];
	}
	return orbit;
}

function computeInnerOrbit(orbit, p, q, w, h, t, n) {
	let orbitLength = orbit.length;
	let innerOrbit = [];
	let iterate = [];

	if (orbitLength>1) {
		for (let i=0; i<orbitLength-1; i++) {
			innerOrbit[i] = ponceletInnerIterate(orbit[i][0], orbit[i][1], orbit[i+1][0], orbit[i+1][1], p, q, w, h, t, n);
		}
	}
	return innerOrbit;
}

function drawOrbit(orbit) {
	let orbitLength = orbit.length;
	push();
		stroke(0,0,0);
		fill(0,0,0);
		for (let i=0; i<orbitLength; i++) {
			strokeWeight(sW2);
			if (i<orbitLength-1) {
				line(orbit[i][0],orbit[i][1],orbit[i+1][0],orbit[i+1][1]);
			}
			strokeWeight(sW1);		
			circle(orbit[i][0],orbit[i][1],0.02);
		}
	pop();
}

function drawInnerOrbit(orbit) {
	let orbitLength = orbit.length;
	push();
		stroke(255,0,0);
		fill(255,0,0);
		for (let i=0; i<orbitLength; i++) {
			strokeWeight(sW1);		
			circle(orbit[i][0],orbit[i][1],0.02);
		}
	pop();
}

function drawCicumscribingSquare() {
	push();
		stroke(200,200,200);
		strokeWeight(sW1);
		noFill();
		beginShape();
			vertex(1,-1);
			vertex(1,1);
			vertex(-1,1);
			vertex(-1,-1);
		endShape(CLOSE);	
	pop();
}

function drawSquareIntersectionsLine(intersections) {
	let intersectionsLength = intersections.length;
	push();
		stroke(0,0,0);
		fill(0,0,0);
		strokeWeight(sW2);
		for (let i=0; i<intersectionsLength-1; i++) {
			line(intersections[i][0],intersections[i][1],intersections[i][2],intersections[i][3]);	
			circle(intersections[i][0],intersections[i][1],0.02);
			circle(intersections[i][2],intersections[i][3],0.02);
		}
	pop();
}

function drawSquareIntersectionsCurve(intersections) {
	let intersectionsLength = intersections.length;
	push();
		stroke(0,0,0);
		fill(0,0,0);
		strokeWeight(sW2);
		beginShape();
			for (let i=0; i<intersectionsLength-1; i++) {				
				curveVertex(intersections[i][0],intersections[i][1]);
				curveVertex(intersections[intersectionsLength-i-2][0],intersections[intersectionsLength-i-2][1]);
				circle(intersections[i][0],intersections[i][1],0.02);
				circle(intersections[i][2],intersections[i][3],0.02);
			}
		endShape();
	pop();
}

function drawSquareIntersectionsBezier(intersections) {
	let intersectionsLength = intersections.length;
	if (intersectionsLength>=4) {
		push();
			stroke(0,0,0);
			strokeWeight(sW2);
			beginShape();
				for (let i=0; i<intersectionsLength-3; i+=4) {
					noFill();
					bezier(intersections[i][0],intersections[i][1],intersections[i+2][0],intersections[i+2][1],intersections[i+1][0],intersections[i+1][1],intersections[i+3][0],intersections[i+3][1]);
					fill(0,0,0);
					circle(intersections[i][0],intersections[i][1],0.02);
					circle(intersections[i][2],intersections[i][3],0.02);
				}
			endShape();
		pop();
	}
}

function drawBezierFlower1(orbit, innerOrbit) {
	let orbitLength = orbit.length;
	let innerOrbitLength = innerOrbit.length;
	let cp0, cp1, cp2, cp3;
	let tgntControlPoints = [];

	if (innerOrbitLength>0) {
		push();
			stroke(0,0,0);
			strokeWeight(sW2);
			beginShape();
				for (let i=0; i<innerOrbitLength; i++) {
					cp0 = orbit[i];
					cp3 = orbit[i+1];
					if (i%2==0) {
						tgntControlPoints = tangentControlPoints(cp0[0],cp0[1],0.5);
						cp1 = tgntControlPoints[1];
						cp2 = innerOrbit[i];
					} else {
						tgntControlPoints = tangentControlPoints(cp3[0],cp3[1],0.5);
						cp1 = innerOrbit[i];
						cp2 = tgntControlPoints[0];
					}
										
					noFill();
					bezier(cp0[0],cp0[1],cp1[0],cp1[1],cp2[0],cp2[1],cp3[0],cp3[1]);
				}
			endShape(); 
		pop();
	}
}

function drawBezierFlower2(orbit, innerOrbit) {
	let orbitLength = orbit.length;
	let innerOrbitLength = innerOrbit.length;
	let cp0, cp1, cp2, cp3;
	// let tgntControlPoints = [];
	let tgntControlPoints0 = [];
	let tgntControlPoints1 = [];

	if (innerOrbitLength>0) {
		push();
			stroke(0,0,0);
			strokeWeight(sW2);
			beginShape();
				// fill(0,0,255,100);
				noFill();
				for (let i=0; i<innerOrbitLength; i++) {
					cp0 = orbit[i];
					cp3 = orbit[i+1];
					tgntControlPoints0 = tangentControlPoints(cp0[0],cp0[1],1.0);
					tgntControlPoints1 = tangentControlPoints(cp3[0],cp3[1],1.0);
					cp1 = tgntControlPoints0[1];
					cp2 = tgntControlPoints1[0];					
										
					bezier(cp0[0],cp0[1],cp1[0],cp1[1],cp2[0],cp2[1],cp3[0],cp3[1]);
				}
				// cp0 = orbit[orbitLength-1];
				// cp3 = orbit[0];
				// tgntControlPoints0 = tangentControlPoints(cp0[0],cp0[1],0.25);
				// tgntControlPoints1 = tangentControlPoints(cp3[0],cp3[1],0.25);
				// cp1 = tgntControlPoints0[1];
				// cp2 = tgntControlPoints1[0];		
				// bezier(cp0[0],cp0[1],cp1[0],cp1[1],cp2[0],cp2[1],cp3[0],cp3[1]);
			endShape(); 
		pop();
	}
}

function drawBezierFlower3(orbit, innerOrbit) {
	let orbitLength = orbit.length;
	let innerOrbitLength = innerOrbit.length;
	let p0, p1;
	let cp0, cp1, cp2, cp3, cp4, cp5, cp6, cp7;
	let tgntControlPointsOuterPre = [];
	let tgntControlPointsOuterPost = [];
	let tgntControlPointsInner = [];

	if (innerOrbitLength>0) {
		push();
			stroke(0,0,0);
			strokeWeight(sW2);
			beginShape();
				for (let i=0; i<innerOrbitLength; i++) {
					p0 = orbit[i];
					p1 = orbit[i+1];
					cp0 = p0;
					cp7 = p1;
					cp3 = innerOrbit[i];
					cp4 = innerOrbit[i];
					tgntControlPointsOuterPre = tangentControlPoints(p0[0],p0[1],0.5);
					tgntControlPointsOuterPost = tangentControlPoints(p1[0],p1[1],0.5);
					tgntControlPointsInner = tangentControlPointsInner(innerOrbit[i][0],innerOrbit[i][1],p0[0],p0[1],p1[0],p1[1],0.25);
					cp1 = tgntControlPointsOuterPre[1];
					cp6 = tgntControlPointsOuterPost[0];
					cp2 = tgntControlPointsInner[0];
					cp5 = tgntControlPointsInner[1];				
										
					noFill();
					bezier(cp0[0],cp0[1],cp1[0],cp1[1],cp2[0],cp2[1],cp3[0],cp3[1]);
					bezier(cp4[0],cp4[1],cp5[0],cp5[1],cp6[0],cp6[1],cp7[0],cp7[1]);								
				}
			endShape(); 
		pop();
	}
}

function drawBezierFlower4(orbit, innerOrbit) {
	let orbitLength = orbit.length;
	let innerOrbitLength = innerOrbit.length;
	let p0, p1;
	let cp0, cp1, cp2, cp3, cp4, cp5, cp6, cp7;
	let tgntControlPointsOuterPre = [];
	let tgntControlPointsOuterPost = [];
	let tgntControlPointsInner = [];

	if (innerOrbitLength>0) {
		push();
			stroke(0,0,0);
			strokeWeight(sW2);
			beginShape();
				for (let i=0; i<innerOrbitLength; i++) {
					p0 = orbit[i];
					p1 = orbit[i+1];
					cp0 = p0;
					cp7 = p1;
					cp3 = innerOrbit[i];
					cp4 = innerOrbit[i];
					tgntControlPointsOuterPre = tangentControlPoints(p0[0],p0[1],1.0);
					tgntControlPointsOuterPost = tangentControlPoints(p1[0],p1[1],1.0);
					tgntControlPointsInner = tangentControlPointsInner(innerOrbit[i][0],innerOrbit[i][1],p0[0],p0[1],p1[0],p1[1],0.5);
					cp1 = tgntControlPointsOuterPre[1]; // Note this (Pre and Post swapped)!
					cp6 = tgntControlPointsOuterPost[0]; // Note this (Pre and Post swapped)!
					cp2 = tgntControlPointsInner[1];
					cp5 = tgntControlPointsInner[0];				
										
					noFill();
					bezier(cp0[0],cp0[1],cp1[0],cp1[1],cp2[0],cp2[1],cp3[0],cp3[1]);
					bezier(cp4[0],cp4[1],cp5[0],cp5[1],cp6[0],cp6[1],cp7[0],cp7[1]);
				}
			endShape(); 
		pop();
	}
}

function drawBezierFlower5(orbit, innerOrbit) {
	let orbitLength = orbit.length;
	let innerOrbitLength = innerOrbit.length;
	let p0, p1;
	let cp0, cp1, cp2, cp3, cp4, cp5, cp6, cp7;
	let tgntControlPointsOuterPre = [];
	let tgntControlPointsOuterPost = [];
	let tgntControlPointsInner = [];

	if (innerOrbitLength>0) {
		push();
			stroke(0,0,0);
			strokeWeight(sW2);
			beginShape();
				for (let i=0; i<innerOrbitLength; i++) {
					p0 = orbit[i];
					p1 = orbit[i+1];
						cp3 = innerOrbit[i];
						cp4 = innerOrbit[i];					
						tgntControlPointsOuterPre = tangentControlPoints(p0[0],p0[1],0.25);
						tgntControlPointsOuterPost = tangentControlPoints(p1[0],p1[1],0.25);
						tgntControlPointsInner = tangentControlPointsInner(innerOrbit[i][0],innerOrbit[i][1],p0[0],p0[1],p1[0],p1[1],0.25);
						cp2 = tgntControlPointsInner[0];
						cp5 = tgntControlPointsInner[1];
					if (i%2==0) {
						cp0 = p0;
						cp6 = p1;
						cp1 = tgntControlPointsOuterPre[1];
						cp7 = tgntControlPointsOuterPost[0];
					} else {
						cp0 = tgntControlPointsOuterPre[0];
						cp1 = p0;
						cp6 = tgntControlPointsOuterPost[1];
						cp7 = p1;
					}
					
					
					
									
										
					noFill();
					bezier(cp0[0],cp0[1],cp1[0],cp1[1],cp2[0],cp2[1],cp3[0],cp3[1]);
					bezier(cp4[0],cp4[1],cp5[0],cp5[1],cp6[0],cp6[1],cp7[0],cp7[1]);
				}
			endShape(); 
		pop();
	}
}

function computeAngularDurations(orbit) {
	let orbitLength = orbit.length;
	let durations = [];
	for (let i=0; i<orbitLength; i++) {
		if (i<orbitLength-1) {
			durations[i] = angularDistance(orbit[i][0],orbit[i][1],orbit[i+1][0],orbit[i+1][1]);
		} else {
			durations[i] = angularDistance(orbit[0][0],orbit[0][1],orbit[orbitLength-1][0],orbit[orbitLength-1][1]);
		}		
	}
	return durations;	
}

function computeEuclideanDurations(orbit) {
	let orbitLength = orbit.length;
	let durations = [];
	for (let i=0; i<orbitLength; i++) {
		if (i<orbitLength-1) {
			durations[i] = euclideanDistance(orbit[i][0],orbit[i][1],orbit[i+1][0],orbit[i+1][1]);
		} else {
			durations[i] = euclideanDistance(orbit[0][0],orbit[0][1],orbit[orbitLength-1][0],orbit[orbitLength-1][1]);
		}	
	}
	return durations;	
}

function angularDistance(x0,y0,x1,y1) {
	return abs(atan2(y0,x0) - atan2(y1,x1));
}

function euclideanDistance(x0,y0,x1,y1) {
	return sqrt(pow(x1-x0,2) + pow(y1-y0,2));
}

function computeSquareIntersections(orbit) {
	let orbitLength = orbit.length;
	let intersections = [];
	let intersectionPoints = [];
	for (let i=0; i<orbitLength; i++) {
		if (i<orbitLength-1) {
			intersectionPoints = squareIntersectionPoints(orbit[i][0],orbit[i][1],orbit[i+1][0],orbit[i+1][1]);
			intersections[i] = [];
			intersections[i][0] = intersectionPoints[0][0];
			intersections[i][1] = intersectionPoints[0][1];
			intersections[i][2] = intersectionPoints[1][0];
			intersections[i][3] = intersectionPoints[1][1];
		} else {
			intersectionPoints = squareIntersectionPoints(orbit[0][0],orbit[0][1],orbit[orbitLength-1][0],orbit[orbitLength-1][1]);
			intersections[i] = [];
			intersections[i][0] = intersectionPoints[0][0];
			intersections[i][1] = intersectionPoints[0][1];
			intersections[i][2] = intersectionPoints[1][0];
			intersections[i][3] = intersectionPoints[1][1];			
		}	
	}
	return intersections;
}

function squareIntersectionPoints(x0,y0,x1,y1) {
	let ts = [];
	let tsSorted = [];
	let ps = [];
	let intersectionPoints = [];
	intersectionPoints[0] = [];
	intersectionPoints[1] = [];
	let leastNegt = -Infinity;
	let leastPost = Infinity;
	let leastNegtIndex = -1;
	let leastPostIndex = -1;	

	ts[0] = -((-2 + x0 + x1)/(x0 - x1));
	ps[0] = [1,(y0 - x1*y0 + (-1 + x0)*y1)/(x0 - x1)];
	ts[1] = -((-2 + y0 + y1)/(y0 - y1));
	ps[1] = [(x0 + x1*(-1 + y0) - x0*y1)/(y0 - y1),1];
	ts[2] = -((2 + x0 + x1)/(x0 - x1));
	ps[2] = [-1,(-((1 + x1)*y0) + (1 + x0)*y1)/(x0 - x1)];
	ts[3] = -((2 + y0 + y1)/(y0 - y1));
	ps[3] = [(x1*(1 + y0) - x0*(1 + y1))/(y0 - y1),-1];

	for (let i=0; i<4; i++) {
		if (abs(ts[i]) != Infinity && ts[i] != NaN) { 
			if (ts[i]<=0) {
				if (ts[i]>leastNegt) {
					leastNegt = ts[i];
					leastNegtIndex = i;
				}
			} else {
				if (ts[i]<leastPost) {
					leastPost = ts[i];
					leastPostIndex = i;
				}				
			}
		}
	}

	intersectionPoints[0][0] = ps[leastPostIndex][0];
	intersectionPoints[0][1] = ps[leastPostIndex][1];
	intersectionPoints[1][0] = ps[leastNegtIndex][0];
	intersectionPoints[1][1] = ps[leastNegtIndex][1];

	return intersectionPoints;

}

function transformCoordinatesWorld(x, y) {
	let transformedCoordinates = [];
	transformedCoordinates[0] = (x-width/2)/scaleX;
	transformedCoordinates[1] = (-y+height/2)/scaleY;	
	return transformedCoordinates;	
}

function transformCoordinatesEllipse(x, y, p, q, t) {
	let transformedCoordinates = [];
	let transformedCoordinatesWorld = transformCoordinatesWorld(x,y);
	let xWorld = transformedCoordinatesWorld[0];
	let yWorld = transformedCoordinatesWorld[1];

	transformedCoordinates[0] = cos(-t)*(xWorld-p)-sin(-t)*(yWorld-q);
	transformedCoordinates[1] = sin(-t)*(xWorld-p)+cos(-t)*(yWorld-q);	
	return transformedCoordinates;
}

function ellipseInsideCircle(p, q, w, h) {
	let longerAxisLength = max(w, h);
	let distanceCentreOfEllipseToCircle = abs(sqrt(pow(p,2)+pow(q,2))-1);
	if (longerAxisLength<distanceCentreOfEllipseToCircle) {	
		return true;
	} else {
		return false;
	}
}

// Returns the clockwise and anticlockwise tangential control
// points - in that order - for iterates on the unit circle
function tangentControlPoints(x, y, s) {
	let controlPoint = [];
	let controlPoints = [];

	controlPoints[0] = [];
	controlPoints[0][0] = x+s*y/sqrt(pow(x,2)+pow(y,2));
	controlPoints[0][1] = y-s*x/sqrt(pow(x,2)+pow(y,2));

    controlPoints[1] = [];
	controlPoints[1][0] = x-s*y/sqrt(pow(x,2)+pow(y,2));
	controlPoints[1][1] = y+s*x/sqrt(pow(x,2)+pow(y,2));

	return controlPoints;
}

// Returns the tangential control points - in that order - for 
// iterates on the inner ellipse.
// These can be 'back' and 'font' depening on the choice of
// circle iterate (c0,c1).
function tangentControlPointsInner(e0, e1, c00, c01, c10, c11, s) {
	let controlPoint = [];
	let controlPoints = [];

	// 'back'
	controlPoints[0] = [];
	controlPoints[0][0] = e0+s*(c00-e0)/sqrt(pow(c00-e0,2)+pow(c01-e1,2));
	controlPoints[0][1] = e1+s*(c01-e1)/sqrt(pow(c00-e0,2)+pow(c01-e1,2));

	// 'front'
	controlPoints[1] = [];
	controlPoints[1][0] = e0+s*(c10-e0)/sqrt(pow(c10-e0,2)+pow(c11-e1,2));
	controlPoints[1][1] = e1+s*(c11-e1)/sqrt(pow(c10-e0,2)+pow(c11-e1,2));

	return controlPoints;
}

export function mousePressed() {
	movableEllipse.clicked(mouseX,mouseY);
}

export function mouseReleased() {
	movableEllipse.stopDragging();
}

export function keyPressed() {
	if (key == 'l') {
		linesOn = !linesOn;		
	}

	if (key == 'c') {
		curvesOn = !curvesOn;
	}

	if (key == 'b') {
		beziersOn = !beziersOn;
	}

	if (key == '0') {
		bezier0On = !bezier0On;
	}

	if (key == '1') {
		bezier1On = !bezier1On;
	}	

	if (key == '2') {
		bezier2On = !bezier2On;
	}	

	if (key == '3') {
		bezier3On = !bezier3On;
	}	

	if (key == '4') {
		bezier4On = !bezier4On;
	}	

	if (key == '5') {
		bezier5On = !bezier5On;
	}						

	if (key == 'u') {
		unitCircleOn = !unitCircleOn;
	}

	if (key == 'e') {
		ellipseOn = !ellipseOn;
	}

	if (key == 'p') {
		ponceletOrbitOn = !ponceletOrbitOn;
	}

	if (key == 's') {
		squareOn = !squareOn;
	}

}
