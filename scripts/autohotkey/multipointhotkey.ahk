#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
; #Warn  ; Enable warnings to assist with detecting common errors.
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.

global cxpos := 318, cypos := 215, count := 1

/*
  Wait for a window to be created, returns 0 on timeout and ahk_id otherwise
  Parameter are the same as WinWait, see http://ahkscript.org/docs/commands/WinWait.htm
  Forum: http://ahkscript.org/boards/viewtopic.php?f=6&t=1274&p=8517#p8517
*/
WinWaitCreated( WinTitle:="", WinText:="", Seconds:=0, ExcludeTitle:="", ExcludeText:="" ) {
    ; HotKeyIt - http://ahkscript.org/boards/viewtopic.php?t=1274
    static Found := 0, _WinTitle, _WinText, _ExcludeTitle, _ExcludeText 
         , init := DllCall( "RegisterShellHookWindow", "UInt",A_ScriptHwnd )
         , MsgNum := DllCall( "RegisterWindowMessage", "Str","SHELLHOOK" )
         , cleanup:={base:{__Delete:"WinWaitCreated"}}
  If IsObject(WinTitle)   ; cleanup
    return DllCall("DeregisterShellHookWindow","PTR",A_ScriptHwnd)
  else if (Seconds <> MsgNum){ ; User called the function
    Start := A_TickCount, _WinTitle := WinTitle, _WinText := WinText
    ,_ExcludeTitle := ExcludeTitle, _ExcludeText := ExcludeText
    ,OnMessage( MsgNum, A_ThisFunc ),  Found := 0
    While ( !Found && ( !Seconds || Seconds * 1000 < A_TickCount - Start ) ) 
      Sleep 16                                                         
    Return Found,OnMessage( MsgNum, "" )
  }
  If ( WinTitle = 1   ; window created, check if it is our window
    && ExcludeTitle = A_ScriptHwnd
    && WinExist( _WinTitle " ahk_id " WinText,_WinText,_ExcludeTitle,_ExcludeText))
    WinWait % "ahk_id " Found := WinText ; wait for window to be shown
}



Numpad4::
CoordMode, Mouse, Screen
MouseGetPos, xpos, ypos ;grab the current mouse position
MsgBox, x %xpos% y %ypos% ;output the mouse position to a message box (note % signs around variables)
return

y::
count := (count+1)
Tooltip, %count%
return


/* ; this is a comment block

Numpad2::
CoordMode, Pixel, Screen
CoordMode, Mouse, Screen
MouseGetPos, xpos, ypos
PixelGetColor, pix, 799, 90
MsgBox, %pix%
return

*/ ; end of comment block


Numpad2:: ;navigate between clusters
if (count == 22){
cxpos := 318
cypos := 341
count := 1
}

CoordMode, Mouse, Screen
MouseMove, cxpos, cypos

cxpos := (cxpos+104)
count := (count+1)
Sleep, 500
Click Left
Click Left
return

Numpad6:: ;get out of "no relevant markers in heap to add" dialogue box
CoordMode, Mouse, Screen
MouseMove, 1470, 783
Click Left
MouseMove, 2522, 39
Click Left
return



Numpad7::
CoordMode, Mouse, Screen
CoordMode, Pixel, Screen
MouseMove, 36, 108
Click Left ;click order markers button
PixelGetColor, pix, 312, 172

PixelSearch, xpix, ypix, 338, 170, 523, 247, 0xE1FFFF
while xpix == ""
{
	;PixelGetColor, pix, 343, 193
	PixelSearch, xpix, ypix, 338, 170, 523, 247, 0xE1FFFF
	Sleep, 10
}
Sleep, 100

/* ;skipping this code for now as not controlling for monotony

MouseMove, 736, 1241 ;move to monotonic button
Sleep, 200
Click Left
Sleep, 200
MouseMove, 1321, 795 ;move to yes button on codominant marker prompt
Click Left
;PixelGetColor, pix, 918, 1235
;while pix != 0x8A3800
;{
;	PixelGetColor, pix, 918, 1235
;	Sleep, 10
;}
Sleep, 1000
while A_Cursor == "Wait"
{
	Sleep, 10
}

*/

Sleep, 300
MouseMove, 470, 98
Click Left
MouseMove, 746, 94
Send !{f4} ; send alt + f4 keypress
Sleep, 500
;insert additional markers and close cluster window
MouseMove, 344, 77
Click Left
Sleep, 200
MouseMove, 394, 101
Click Left
Sleep, 200
MouseMove, 892, 710
Click Left
MouseMove, 2519, 37
Sleep, 200
Click Left
MouseMove, 1199, 284



;Sleep, 500
;PixelSearch, xpix2, ypix2, 428, 96, 1025, 105, 0x7D89E3
;MouseMove, xpix2, ypix2


;MouseGetPos, xpos, ypos
;WinGetPos, x, y, width, height, A
;MsgBox, mouse x %xpos% mouse y %ypos% win wid %width% win height %height%
;MouseMove, x, y





;Sleep, 2000




;WinWaitCreated()
;MsgBox, a new window has opened
;Return
return

Numpad5::
WinGetPos, x, y, width, height, A
MsgBox, width %width% height %height%
return

Numpad1::
CoordMode, Mouse, Screen
MouseMove, 344, 77
Click Left
Sleep, 200
MouseMove, 394, 101
Click Left
Sleep, 200
MouseMove, 892, 710
Click Left
MouseMove, 2519, 37
Sleep, 200
Click Left
MouseMove, 1199, 284
return

Numpad9::
Tooltip, %A_Cursor%
return


^s::
CoordMode, Mouse, Screen
MouseGetPos, xpos, ypos
MouseMove, 163, 45
Click Left
MouseMove, xpos, ypos
return

