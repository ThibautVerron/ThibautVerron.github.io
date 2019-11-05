function showMenuByWidth() {
    var x = document.getElementById("menu_dropdown");
    if (x.className.indexOf("w3-show") == -1) {
        x.className += " w3-show";
    } else { 
        x.className = x.className.replace(" w3-show", "");
    }
}

// window.onload = function () {
//     document.getElementById("get-apropos-horiz").onclick = function () {
//         console.log("toto")
//         document.getElementById("apropos").scrollIntoView(true);
//     };
//     // document.getElementById("get-apropos-horiz").color = red;
// }

// https://www.the-art-of-web.com/javascript/remove-anchor-links/
window.addEventListener("DOMContentLoaded", function(e) {
  var links = document.getElementsByTagName("A");
  for(var i=0; i < links.length; i++) {
    if(!links[i].hash) continue;
    if(links[i].origin + links[i].pathname != self.location.href) continue;
    (function(anchorPoint) {
      links[i].addEventListener("click", function(e) {
        anchorPoint.scrollIntoView(true);
        e.preventDefault();
      }, false);
    })(document.getElementById(links[i].hash.replace(/#/, "")));
  }
}, false);
