export const INTERVALS = {
	major: [0, 2, 4, 5, 7, 9, 11],
	minor: [0, 2, 3, 5, 7, 8, 10]
}

export const CHROMATIC_SCALE = ['c', 'c#', 'd', 'd#', 'e' , 'f', 'f#', 'g', 'g#', 'a', 'a#', 'b']


export const SCALE = {
	major: buildScaleLookup('major'),
	minor: buildScaleLookup('minor'),
	custom: intervals => buildScaleLookup(intervals)
}

export function buildScaleLookup(type) {
	return Object.fromEntries(CHROMATIC_SCALE.map(n => [n, buildScale(n, type)]));
}

export function buildScale(root, type = []) {
	const notes = sortByRootNote(root);
	if(Array.isArray(type)){
		return [...type].slice(0, 7).map(i => notes[i])
	}
	return INTERVALS[type].map(i => notes[i])
}

export function noteValue(note, octave) {
	return octave * 12 + CHROMATIC_SCALE.indexOf(note)
}

export function noteString(value) {
	return `${CHROMATIC_SCALE[value % 12]}${Math.floor(value/12)}`
}

export function buildChord(note, octave, type, count) {
	const value = noteValue(note, octave)
	if(count <= 4) {
		return chordNotes(value, count, [...INTERVALS[type]])
	}
	const [ ...intervals ] = growIntervals(INTERVALS[type], Math.ceil(count/4))
	return chordNotes(value, count, intervals)
}

export function chordNotes(noteValue, noteCount, intervals) {
	return intervals
		.filter((_, i) => i % 2 === 0)
		.slice(0, noteCount)
		.map(i => midiToFreq(noteValue+i))
}

function growIntervals(initial, noteCount) {
	let count = 1;
	const intervals = [initial]
	while(count < noteCount) {
		intervals.push(initial.map(i => 12 * count + i ))
		count++
	}
	return intervals.flat() 
}

function sortByRootNote(root) {
	if(root === 'c') return CHROMATIC_SCALE
	const index = CHROMATIC_SCALE.indexOf(root)
	const left = CHROMATIC_SCALE.slice(0, index)
	const right = CHROMATIC_SCALE.slice(index)
	return [ ...right, ...left ];
}
