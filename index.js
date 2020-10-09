import * as sketch from './sketch.js';

for (let k in sketch) {
	window[k] = sketch[k];
}
