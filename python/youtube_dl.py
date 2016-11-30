import youtube_dl
import os 

options = {'outtmpl': '%(id)s'} 

options = {
    'format': 'bestaudio/best', # choice of quality
    'extractaudio' : True,      # only keep the audio    # convert to mp3 
    'outtmpl': '%(id)s',  
    'audioformat' : "mp3",      # name the file the ID of the video
    'noplaylist' : True,        # only download single song, not playlist
}



with youtube_dl.YoutubeDL(options) as ydl:
    ydl.download(['https://www.youtube.com/watch?v=GXoZLPSw8U8'])

