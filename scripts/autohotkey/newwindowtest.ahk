#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
; #Warn  ; Enable warnings to assist with detecting common errors.
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.

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

;Numpad1::
;CoordMode, Mouse, Screen
;MouseGetPos, xpos, ypos ;grab the current mouse position
;MsgBox, x %xpos% y %ypos% ;output the mouse position to a message box (note % signs around variables)
;return

;Numpad2::
;CoordMode, Pixel, Screen
;PixelGetColor, pix, 343, 193
;MsgBox, %pix%
;return

Numpad0::
CoordMode, Mouse, Screen
CoordMode, Pixel, Screen
MouseMove, 36, 108
Click Left
PixelGetColor, pix, 312, 172
while pix != 0xE1FFFF
{
	PixelGetColor, pix, 343, 193
	Sleep, 10
}
Sleep, 100

MouseMove, 736, 1241
Sleep, 200
Click Left
Sleep, 200
MouseMove, 1334, 793
Click Left

;WinWaitCreated()
;MsgBox, a new window has opened
;Return
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
return




