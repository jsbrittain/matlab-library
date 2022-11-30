<!--- Hide from tired old browsers

function floater(num)
{ var floater = window.open(num+"floater.html", "BackgroundMusic", "width=140,height=30,Scrollbars");
}

function jumpBox(list,num) 
 {
  if(list.elements[0].options[0].value==666)
    {list.elements[0].length=19;
    list.elements[0].options[0].selected=1;
    list.elements[0].options[0].text='[  WarpMe  ]';
    list.elements[0].options[0].value='BLANK';
    list.elements[0].options[1].text='Main';
    list.elements[0].options[1].value=num+'default.html';
    list.elements[0].options[2].text='Personal';
    list.elements[0].options[2].value=num+'personal/personal.html';
    list.elements[0].options[3].text='Tutorials';
    list.elements[0].options[3].value=num;
    list.elements[0].options[4].text='Projects';
    list.elements[0].options[4].value=num+'projects/projects.html';
    list.elements[0].options[5].text='Links';
    list.elements[0].options[5].value=num+'links/links.html';
    list.elements[0].options[6].text='WebRings';
    list.elements[0].options[6].value=num+'webrings/webring.html'; 
    list.elements[0].options[7].text='Files';
    list.elements[0].options[7].value=num+'files/filelist.html';
    list.elements[0].options[8].text='Reviews';
    list.elements[0].options[8].value=num+'bookreviews/bookreviews.html';
    list.elements[0].options[9].text='Banners';
    list.elements[0].options[9].value=num+'banners/banners.html';
    list.elements[0].options[10].text='News';
    list.elements[0].options[10].value=num+'news/News.html';
    list.elements[0].options[11].text='Feedback';
    list.elements[0].options[11].value=num+'feedback/Feedback.html';
    list.elements[0].options[12].text='Copyrights';
    list.elements[0].options[12].value=num+'copyrights/copyrights.html';
    list.elements[0].options[13].text='Email';
    list.elements[0].options[13].value='mailto:deltener@mindtremors.com';
    list.elements[0].options[14].text='Subscribe';
    list.elements[0].options[14].value=num+'mailing%20list/mailinglist.html';
    list.elements[0].options[15].text='Help';
    list.elements[0].options[15].value=num+'help/help.html';
    list.elements[0].options[16].text='Message Board';
    list.elements[0].options[16].value=num+'cgi-bin/messageboard.pl?Mode=DisplaySubject';
    list.elements[0].options[17].text='Beta Testing ';  //'Search'; 
    list.elements[0].options[17].value=num+'betatesting/betatesting.html';
    list.elements[0].options[18].text='My Mission';
    list.elements[0].options[18].value=num+'mission/mission.html';
   }
  else if(list.elements[0].options[list.elements[0].selectedIndex].value =='ReLoad')
   { list.elements[0].options[0].value=666;
     list.elements[0].options[1].selected=1;
     jumpBox(list,num);
   }
  else if(list.elements[0].options[list.elements[0].selectedIndex].text =='Tutorials')
   {  list.elements[0].length=9;
      list.elements[0].options[0].selected=1;
      list.elements[0].options[0].text='Pick One';
      list.elements[0].options[0].value='BLANK';
      list.elements[0].options[1].text='<--- Back ';
      list.elements[0].options[1].value='ReLoad';
      list.elements[0].options[2].text=' ';
      list.elements[0].options[2].value='BLANK';
      list.elements[0].options[3].text='Sound';
      list.elements[0].options[3].value=num+'tutorials/sound%20programming/soundprogramming.html';
      list.elements[0].options[4].text='Graphics';
      list.elements[0].options[4].value=num+'tutorials/graphics%20programming/graphicsprogramming.html';
      list.elements[0].options[5].text='Pmode';
      list.elements[0].options[5].value=num+'tutorials/djgpp%20programming/pmode.html';
      list.elements[0].options[6].text='Interrupt';
      list.elements[0].options[6].value=num+'tutorials/interrupt%20programming/interruptprogramming.html';
      list.elements[0].options[7].text='Input Devices';
      list.elements[0].options[7].value=num+'tutorials/input%20device%20programming/inputdevices.html';
      list.elements[0].options[8].text='C++';
      list.elements[0].options[8].value=num+'tutorials/C++/C++.html';
   }
  else 
   { if((list.elements[0].options[list.elements[0].selectedIndex].value) != 'BLANK')
        { location.href = list.elements[0].options[list.elements[0].selectedIndex].value;
          list.elements[0].selectedIndex=0;
        }
   } 
 }

 // end hiding --->