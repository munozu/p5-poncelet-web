import * as audio from './audio.js'
import * as poncelet from './poncelet.js'

let audioData;
let ponceletData;

export function setup() {
	poncelet.setup(audioData, setPonceletData)
	audio.setup(ponceletData, setAudioData)
}

export function draw() {
	poncelet.draw(audioData, setPonceletData)
	audio.draw(ponceletData, setAudioData)
	// if(frameCount % 100 === 0) {
	// 	console.log(audioData)
	// }
}

function setAudioData(new_data) {
	audioData = new_data;
}

function setPonceletData(new_data) {
	ponceletData = new_data;
}

export { mousePressed, mouseReleased, keyPressed } from './poncelet.js'
