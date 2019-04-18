
// ------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------- //
// 										GUITAR EFFECTS ON BELA										 //
// 								Assignment 3 - Real-Time DSP (ECS732P)							     //
//									Kranthi Yanamandra (170797647)									 //
//																									 //
// ------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------- //

#include <Bela.h>
#include <cmath>
#define DELAY_BUFFER_SIZE 44100
#include <WriteFile.h>

// Variables for reading audio input
float input = 0.0;
float input_ = 0.0;


float gInverseSampleRate;

// Output variables for different effects
float output_1 = 0.0;
float output_2 = 0.0;
float output_3 = 0.0;
float output_4 = 0.0;
float output_5 = 0.0;


// Low Frequency Oscillators' variables
float gLfoPhase_vibrato;
float glfoPhase_tremolo;
float glfoPhase_ringmod;
float gLfoFrequency_1 = 0.0;
float gLfoFrequency_2 = 0.0;

// Used for distortion 
float gGain = 0.0;

// Number of audio frames per analog frame
int gAudioFramesPerAnalogFrame;

// Read and Write pointers for delay buffer
float gReadPointer = 0.0;
int gWritePointer = 0;

// Delay buffer
float gDelayBuffer[DELAY_BUFFER_SIZE];


// Variables for delay based effects
float gSweepWidth = 0.01;
float gInterpolatedSample = 0.0;
float gInterpolatedSample_ = 0.0;
float gCurrentDelay = 0.0;
float gDelay = 0.0025;
float gFeedback = 0.5;

// Variables for Ring Modulation
float gCarrierPhase;
float gCarrierFrequency;

bool setup(BelaContext *context, void *userData)
{
	gInverseSampleRate = 1.0 / context->audioSampleRate;
	glfoPhase_ringmod = 0.0;
	glfoPhase_tremolo = 0.0;
	gLfoPhase_vibrato = 0.0;
	
	
	gCarrierPhase = 0.0;
	gCarrierFrequency = 100;

	gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;
	return true;
}

void render(BelaContext *context, void *userData)
{
	for(unsigned int n = 0; n < context->audioFrames; n++)
	{
			
		for(unsigned int channel = 0; channel < context->audioInChannels; channel++)
		{
			// Reading in the audio
			input = audioRead(context, n, channel);
			input_ = audioRead(context, n, channel);
			
			if(!(n % gAudioFramesPerAnalogFrame)) 
			{
				// On even audio samples:
		
		
				// ------------------------------------------- Vibrato ------------------------------------------------ //
			
				// Create LFO for vibrato
				float lfo_vibrato = sinf(2.0 * M_PI * gLfoPhase_vibrato) * 0.5 + 0.5f;
				
				// Keep track and wrap the phase of the sinewave
				gLfoPhase_vibrato += gLfoFrequency_1 * gInverseSampleRate;
				if(gLfoPhase_vibrato > 1.0)
				gLfoPhase_vibrato -= 1.0;
				
				// Adjust frequency of LFO based on potentiometer reading
				gLfoFrequency_1 = map(analogRead(context, n, 1), 0.0, 1.0, 0.0, 10);
				
				// update current delay
				gCurrentDelay = gSweepWidth * lfo_vibrato;
				
				// calculate the read pointer position in terms of the write pointer. 
				gReadPointer = fmodf((float)gWritePointer - (float)(gCurrentDelay * context->audioSampleRate) + (float)(DELAY_BUFFER_SIZE) - 3.0, (float)DELAY_BUFFER_SIZE);
				
				
				// Cubic interpolation
				
				int sample_1 = (int)floorf(gReadPointer);
				int sample_2 = (sample_1 + 1) % DELAY_BUFFER_SIZE;
				int sample_3 = (sample_2 + 1) % DELAY_BUFFER_SIZE;
				int sample_0 = (sample_1 - 1 + DELAY_BUFFER_SIZE) % DELAY_BUFFER_SIZE;
				
				float fraction = gReadPointer - floorf(gReadPointer);
				float fracsquare = fraction * fraction;
				
				float a0 = -0.5f * gDelayBuffer[sample_0] + 1.5f * gDelayBuffer[sample_1] - 1.5f * gDelayBuffer[sample_2] + 0.5f * gDelayBuffer[sample_3];
				float a1 = gDelayBuffer[sample_0] - 2.5f * gDelayBuffer[sample_1] + 2.0f * gDelayBuffer[sample_2] - 0.5f * gDelayBuffer[sample_3];
				float a2 = -0.5f * gDelayBuffer[sample_0] + 0.5f * gDelayBuffer[sample_2];
                float a3 = gDelayBuffer[sample_1];
                
                // calculate interpolated sample
                gInterpolatedSample = a0 * fraction * fracsquare + a1 * fracsquare + a2 * fraction + a3;

				// write to the delay buffer
				gDelayBuffer[gWritePointer] = input;
				
				// house keeping for the write pointer
				if(++gWritePointer >= DELAY_BUFFER_SIZE)
				gWritePointer = 0;
				
				// output for the vibrato effect
				output_1 = gInterpolatedSample;
				
			
			// ------------------------------------------------ Tremolo ------------------------------------------------  //
			
				// Create LFO for tremolo effect
				float lfo_tremolo = sinf(glfoPhase_tremolo) * 0.5 + 0.5f;
				
				// Adjust frequency of LFO based on potentiometer reading.  
				gLfoFrequency_2 = map(analogRead(context, n, 2), 0.0, 1.0, 0.0, 15.0);
				
				// Keep track and wrap the phase of the sinewave.
				glfoPhase_tremolo +=  2.0 * M_PI * gLfoFrequency_2 * gInverseSampleRate;
				if(glfoPhase_tremolo >= 2.0 * M_PI)
				glfoPhase_tremolo -= 2.0 * M_PI;
				
				// output for the tremolo effect
				output_2 =  input * lfo_tremolo;
			
			
			
			// ------------------------------------------ Ring modulation ----------------------------------------- //
			
				// output for the ring modulation effect
				output_3 = input * sinf(2.0 * M_PI * gCarrierPhase);
				
				// update carrier frequency based on potentiometer reading
				gCarrierFrequency = map(analogRead(context, n, 3), 0.0, 1.0, 50, 500);
				
				// Create LFO for carrier. Update carrier and lfo phases. Keep in range 0 to 1
				float lfo_ringmod = sinf(glfoPhase_ringmod) * 0.5 + 0.5f;
				glfoPhase_ringmod += gLfoFrequency_1 * gInverseSampleRate;
				if(glfoPhase_ringmod >= 1)
					glfoPhase_ringmod -= 1.0;
					
				gCarrierPhase += (gCarrierFrequency + gSweepWidth * lfo_ringmod) * gInverseSampleRate;
            	if(gCarrierPhase >= 1.0)
                gCarrierPhase -= 1.0; 
				
				
			
	
			// --------------------------------------------- Distortions ---------------------------------------- //
			
			/* Soft exponential clipping */
				
				// Calculate input gain
				float gain = powf(10.0f, gGain/20.0f);
				
				// apply gain to the input
				float in = input_ * gain;
				
				// update gain by reading potentiometer values
				gGain = map(analogRead(context, channel, 0), 0.0, 1.0, 0.0, 24.0);
				
				// Calculate exponential clipping
				if(in > 0.0)
				{
					output_4 = 1.0 - expf(-in);
				}
				else
				{
					output_4 = - 1.0 + expf(in);
				}
				
				// output for soft exponential clipping
				output_4 /= 2.0;
				
			
			/* Half Wave Rectification */
				
				// output for half wave rectification
		  		output_5 = 0.5 * (fabs(input) + input);
		  		
		  		// Final audio output, combining all the previous outputs from each effect
				audioWrite(context, n, channel, output_1 + output_2 + output_3 + output_4 + output_5);
				
				//audioWrite(context, n, channel, output_1);
				//audioWrite(context, n, channel, output_2);
				//audioWrite(context, n, channel, output_3);
				//audioWrite(context, n, channel, output_4);
				//audioWrite(context, n, channel, output_5);

			
			// ------------------------------------------------ xxxx -------------------------------------------- //
			
			

			}
			
		
	}
}
}
void cleanup(BelaContext *context, void *userData)
{

}