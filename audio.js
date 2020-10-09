import {CHROMATIC_SCALE, SCALE, noteValue, noteString, buildChord} from './musicUtil.js';

const BPM = 120 

const COF = ['a', 'd', 'g', 'c', 'f', 'a#', 'd#', 'g#', 'c#', 'f#', 'b', 'e']

const COF_LOOKUP = Object.fromEntries(COF.map((n,i) => {
	const radians = i * Math.PI/6;
	const coords = [Math.cos(i*Math.PI/6), Math.sin(i*Math.PI/6)];
	const subdivisions = [...Array(6).keys()].map(j => {
		let min = radians - (j + 1) * Math.PI/72;
		let max = radians +  (j + 1) * Math.PI/72;
		return { min, max }
	}); 
	return [
		n, 
		{ radians, coords, subdivisions }
	]
}
))

const COF_RADIANS = Object.values(COF_LOOKUP).map(n => n.radians)
const COF_COORDS =  Object.values(COF_LOOKUP).map(n => n.coords)
const COF_DIVISIONS =  Object.values(COF_LOOKUP).map(n => n.subdivisions)

let data;
let playing = false
let step = 0
let primary_note, secondary_note, primary_chord;
let fft;

export function setup(ponceletData, sendData) {
	if(ponceletData) data = ponceletData 

	const btn = createButton('play/pause').size(96, 40);
	setupPoly()
	setupNoise()
	setupMelody()
	setupComp()
	btn.position(0,0)
	btn.mousePressed(play);
	userStartAudio()
	fft = new p5.FFT()
}


let comp;
function setupComp() {
	comp = new p5.Compressor()
	comp.attack(0.8)
	comp.release(0.5)
	comp.ratio(4)
	comp.threshold(-30)
	comp.drywet(0.5)
}

export function draw(sendData) {
	sendData({
		fft,
		spectrum: fft.analyze(),
		waveform: fft.waveform()
	})
}

let poly_loop, poly_fltr, poly_chord, poly_phrase, poly_rev, poly_step = 0;
function setupPoly() {
	poly_chord = new p5.PolySynth()
	poly_chord.setADSR(4, 3, 0, 8)
	poly_rev = new p5.Reverb()
	poly_rev.process(poly_chord, 8, 0.5)
	poly_rev.drywet(1)
	poly_fltr = new p5.HighPass()
	poly_fltr.process(poly_rev)
	poly_fltr.freq(210)
	poly_fltr.res(2)
	poly_fltr.drywet(1)
	poly_fltr.connect(comp)

	poly_phrase = new p5.Phrase('poly', playPoly, data.ponceletOrbit)
	poly_loop = new p5.Part(null, 2/1)
	poly_loop.setBPM(BPM)
	poly_loop.addPhrase(poly_phrase)
}

function playPoly(time, [x,y]) {
	const scale = data.pVal < 0 ? 'minor': 'major'
	const chord_size = 
		poly_step % 4 === 0 
		? 3 
		: poly_step % 8 === 0 
		? 2 
		: poly_step % 10 === 0 
		? 5 
		: 4
	const a = getAngle(x,y)
	const n = getRootNote(a)
	const diatonic = tonic(a)
	let d_note = SCALE[scale][n][diatonic]
	primary_note = buildChord(d_note, 4, scale, 1)[0];
	secondary_note = buildChord(n, 4, scale, 1)[0];
	const chord = [primary_note/2, ...buildChord(d_note, 4, scale, chord_size)]
	chord.forEach((n, i) => {
		poly_chord.play(n, 0.01, time + (i * 0.017), 1)
	})
	if(Math.random() < 0.5) poly_chord.play(chord.pop() * 2, 0.02, time += 1/16, 0.2)
	primary_chord = chord
	poly_step = (poly_step+1) % 16
}

let nz_loop, nz_env, nz_flt_env, nz_rytm, nz_phrase, nz_fltr, nz_step = 0;
function setupNoise() {
	nz_rytm = new p5.Noise()
	nz_fltr = new p5.BandPass()
	nz_rytm.disconnect()
	nz_rytm.start()
	nz_rytm.amp(0)

	nz_env = new p5.Env();
	nz_env.setADSR(0.001,0.01);
	nz_env.setRange(0.02, 0)

	nz_flt_env = new p5.Env();

	nz_fltr.process(nz_rytm)
	nz_fltr.connect(comp)
	nz_fltr.drywet(1)

	nz_phrase = new p5.Phrase('nz', playNoise, data.euclideanDurations)
	nz_loop = new p5.Part(null, 1/12)
	nz_loop.setBPM(BPM)
	nz_loop.addPhrase(nz_phrase)
}

function playNoise(time, freq) {
	if(!primary_note) return;
	if(nz_step % 1 === 0) {
		const f = random([ primary_note, primary_note * 2, primary_note * 4 ])
		// nz_flt_env.setADSR(0.001 * (freq * 2), 0.01, 0, 0.2);
		nz_flt_env.setADSR(0.002, 0.002);
		nz_flt_env.setRange(primary_note *2, f)
		nz_fltr.freq(nz_flt_env)
		nz_fltr.res(12)
		nz_rytm.pan(random([-0.5, 0.25, 0.25, 0.5]))
		nz_flt_env.play()
		nz_env.play(nz_rytm)
	}
	nz_step = (nz_step+1) % 16;
}

let melody, mel_loop, mel_phrase, mel_delay, mel_fltr, mel_step = 0;
function setupMelody() {
	melody = new p5.MonoSynth()
	melody.setADSR(0.1, 0.2, 0, 0.5)
	mel_delay = new p5.Delay()
	melody.connect(mel_delay)
	mel_fltr = new p5.HighPass()
	mel_fltr.process(mel_delay)
	mel_fltr.freq(100)
	mel_fltr.res(2)
	mel_fltr.drywet(1)
	mel_fltr.connect(comp)
	// source, delayTime (in seconds), feedback, filter frequency

	mel_phrase = new p5.Phrase('mel', playMelody, data.angularDurations)
	mel_loop = new p5.Part(null, 1/8)
	mel_loop.setBPM(BPM)
	mel_loop.addPhrase(mel_phrase)
}

function playMelody(time) {
	if(!primary_note) return
	mel_delay.delayTime(0.4)
	mel_delay.feedback(0.5)

	if(mel_step % 8 < 4) {
		melody.play(primary_chord[mel_step % 3], 0.01, time, 0.1)
	} else {
		melody.play(primary_chord[mel_step % 3] * 2, 0.01, time, 0.1)
	} 

	mel_step = (mel_step+1) % 16 
}

function play() {
	if(!playing) {
		poly_loop.loop(0);
		nz_loop.loop(0);
		mel_loop.loop(0);
		playing = true
	} else {
		poly_loop.stop(0);
		nz_loop.stop(0);
		mel_loop.stop(0);
		poly_step = 0
		mel_step = 0 
		nz_step = 0 
		playing = false
	}
}

function getNextIndex(index, arr) {
	return arr[index+1] || arr[0]
}


function getData(d){
	data = d
	if(poly_phrase){
		poly_phrase.sequence = d.ponceletOrbit
	}
}

function getAngle(x,y){
	const angle = Math.atan2(y, x)  
	return  angle < 0 ? angle + 2 * Math.PI : angle
}

function getRootNote(angle){
	const idx = COF_RADIANS.findIndex((rad, i, arr) => {
		if((i + 1) >= COF_RADIANS.length) return angle >= rad && angle < 2 * Math.PI
		return angle >= rad && angle < COF_RADIANS[i+1]
	})
	return COF[idx] || COF[0]
};

function tonic(a){
	for(let div of COF_DIVISIONS) {
		const index = div.findIndex((d, i) => d.min <= a &&  a < d.max)
		if(index < 0) continue;
		return index
	}
}
